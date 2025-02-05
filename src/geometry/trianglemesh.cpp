#include "eigrene/geometry/trianglemesh.hpp"

namespace eigrene {

index_t TriangleMesh::n_cells() const { return triangles.size(); }

TriangleIndices TriangleMesh::get_triangle_indices(index_t cell_index) const
{
    return triangles[cell_index];
}

bool TriangleMesh::is_in_cell(Vec2 pos, index_t cell_index) const
{
    TriangleIndices tind = triangles[cell_index];
    Triangle t(vertices[tind.ind1], vertices[tind.ind2], vertices[tind.ind3]);

    return t.contains(pos);
}

index_t TriangleMesh::get_cell_index(Vec2 pos) const
{
    auto candidates = spatialindex.query(pos);

    for (index_t i = 0; i < candidates.size(); i++) {
        index_t cell_index = candidates[i];
        if (is_in_cell(pos, cell_index)) {
            return cell_index;
        }
    }

    return -1;
}

real_t TriangleMesh::get_volume(index_t cell_index) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    // This calculation assumes the triangle is stored counter-clockwise
    return 0.5 * (-p2.y * p3.x + p1.y * (-p2.x + p3.x) + p1.x * (p2.y - p3.y) +
                  p2.x * p3.y);
}

// Return the time until the intersection of particle with position pos and
// velocity vel with an edge defined by points p1 and p2. Note, this function
// does not verify whether the returned time is positive or whether the
// intersection lies between the two points. This will be implicitly done later.
inline real_t find_edge_intersection(const Vec2& p1,
                                     const Vec2& p2,
                                     const Vec2& pos,
                                     const Vec2& vel)
{
    // Construct 2x2 matrix equation to find the intersection
    real_t a00 = p2.x - p1.x;
    real_t a10 = p2.y - p1.y;
    real_t a01 = -vel.x;
    real_t a11 = -vel.y;

    real_t b0 = pos.x - p1.x;
    real_t b1 = pos.y - p1.y;

    real_t determinant = a00 * a11 - a01 * a10;

    if (std::abs(determinant) < 1.0e-8) {
        // matrix is singular, no intersection
        // return a negative number because those will be ignored anyway.
        return -1;
    }

    // position of intersection: new_pos = p1 + s * (p2 - p1)
    // real_t s = (a11 * b0 - a01 * b1) / determinant;

    // if (s < 0 or s > 1) {
    //     // intersection is not on the edge, return -1
    //     return -1;
    // }

    return (a00 * b1 - b0 * a10) / determinant;
}

CellCrossing TriangleMesh::get_next_crossing(ParticleGeometry& p) const
{
    assert(is_in_cell(p.position, p.cell_id));
    TriangleIndices t = triangles[p.cell_id];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    // Time till intersection with edge 0
    real_t t0;
    if(p.boundary_id == 0){
        t0 = -1;
    } else {
        t0 = find_edge_intersection(p1, p2, p.position, p.velocity);
    }

    real_t t1;
    if(p.boundary_id == 1){
        t1 = -1;
    } else {
        t1 = find_edge_intersection(p2, p3, p.position, p.velocity);
    }

    real_t t2;
    if(p.boundary_id == 2){
        t2 = -1;
    } else {
        t2 = find_edge_intersection(p3, p1, p.position, p.velocity);
    }

    // t0, t1, and t2 are now the times until the intersection with the edges
    // the true collision time is the smallest positive value of the three
    real_t t_new;
    index_t id_new;
    index_t surface_id;
    if (t0 > 0 and (t0 <= t1 or t1 < 0) and (t0 <= t2 or t2 < 0)) {
        t_new = t0;
        id_new = neighbors[3 * p.cell_id];
        surface_id = 0;
    }
    else if (t1 > 0 and (t1 <= t0 or t0 < 0) and (t1 <= t2 or t2 < 0)) {
        t_new = t1;
        id_new = neighbors[3 * p.cell_id + 1];
        surface_id = 1;
    }
    else if (t2 > 0 and (t2 <= t0 or t0 < 0) and (t2 <= t1 or t1 < 0)) {
        t_new = t2;
        id_new = neighbors[3 * p.cell_id + 2];
        surface_id = 2;
    }
    else {
        throw std::runtime_error("Failed to find next crossing");
    }

    Vec2 pos_new = p.position + p.velocity * t_new;

    CellCrossing crossing = CellCrossing(pos_new, t_new, surface_id);

    if (id_new < 0) {
        // If the id is negative, it is a boundary,
        // so negate it and subtract 1 (zero based) to get the boundary id
        crossing.boundary_condition_id = -id_new - 1;
        return crossing;
    }

    crossing.next_cell_id = id_new;
    // Find the surface id of the next cell
    for (index_t i = 0; i < 3; i++) {
        if (neighbors[3 * id_new + i] == p.cell_id) {
            crossing.next_surface_id = i;
            break;
        }
    }

    return crossing;
}

void make_neighborlist(const vector<TriangleIndices>& triangles,
                       vector<index_t>& neighbors)
{
    if (neighbors.size() != 3 * triangles.size()) {
        neighbors.resize(3 * triangles.size());
    }

    // For now, this is really stupid, we're just manually going through all the
    // triangles and checking the edges
    // A kd tree would massively speed this up
    for (index_t i = 0; i < triangles.size(); i++) {
        TriangleIndices t1 = triangles[i];

        // Initialize to -1, both to keep track if all neighbors of this
        // triangle have been found, but also to keep track of boundaries.
        neighbors[3 * i] = -1;
        neighbors[3 * i + 1] = -1;
        neighbors[3 * i + 2] = -1;

        for (index_t j = 0; j < triangles.size(); j++) {
            if (i == j) continue;

            TriangleIndices t2 = triangles[j];

            bool has_ind1 =
                t2.ind1 == t1.ind1 || t2.ind2 == t1.ind1 || t2.ind3 == t1.ind1;
            bool has_ind2 =
                t2.ind1 == t1.ind2 || t2.ind2 == t1.ind2 || t2.ind3 == t1.ind2;
            bool has_ind3 =
                t2.ind1 == t1.ind3 || t2.ind2 == t1.ind3 || t2.ind3 == t1.ind3;

            if (has_ind1 and has_ind2) {
                neighbors[3 * i] = j;
            }
            if (has_ind2 and has_ind3) {
                neighbors[3 * i + 1] = j;
            }
            if (has_ind3 and has_ind1) {
                neighbors[3 * i + 2] = j;
            }

            // Early exit if all neighbors are found
            if (neighbors[3 * i] != -1 and neighbors[3 * i + 1] != -1 and
                neighbors[3 * i + 2] != -1) {
                break;
            }
        }
    }
}

void make_neighborlist(const vector<Vec2>& vertices,
                       const vector<TriangleIndices>& triangles,
                       const SpatialIndex& spatialindex,
                       vector<index_t>& neighbors)
{
    if (neighbors.size() != 3 * triangles.size()) {
        neighbors.resize(3 * triangles.size());
    }

    for (index_t i = 0; i < triangles.size(); i++) {
        TriangleIndices t1 = triangles[i];

        // Initialize to -1, both to keep track if all neighbors of this
        // triangle have been found, but also to keep track of boundaries.
        neighbors[3 * i] = -1;
        neighbors[3 * i + 1] = -1;
        neighbors[3 * i + 2] = -1;

        // Find the bounding box of the triangle
        Vec2 p1 = vertices[t1.ind1];
        Vec2 p2 = vertices[t1.ind2];
        Vec2 p3 = vertices[t1.ind3];

        real_t min_x = std::min({p1.x, p2.x, p3.x});
        real_t max_x = std::max({p1.x, p2.x, p3.x});
        real_t min_y = std::min({p1.y, p2.y, p3.y});
        real_t max_y = std::max({p1.y, p2.y, p3.y});

        BoundingBox box({min_x, min_y}, {max_x, max_y});

        vector<index_t> candidates = spatialindex.query(box);

        for(index_t j = 0; j < candidates.size(); j++){
            index_t k = candidates[j];
            if(i == k) continue;
            TriangleIndices t2 = triangles[k];

            bool has_ind1 =
                t2.ind1 == t1.ind1 || t2.ind2 == t1.ind1 || t2.ind3 == t1.ind1;
            bool has_ind2 =
                t2.ind1 == t1.ind2 || t2.ind2 == t1.ind2 || t2.ind3 == t1.ind2;
            bool has_ind3 =
                t2.ind1 == t1.ind3 || t2.ind2 == t1.ind3 || t2.ind3 == t1.ind3;

            if (has_ind1 and has_ind2) {
                neighbors[3 * i] = k;
            }
            if (has_ind2 and has_ind3) {
                neighbors[3 * i + 1] = k;
            }
            if (has_ind3 and has_ind1) {
                neighbors[3 * i + 2] = k;
            }

            // Early exit if all neighbors are found
            if (neighbors[3 * i] != -1 and neighbors[3 * i + 1] != -1 and
                neighbors[3 * i + 2] != -1) {
                break;
            }
        }
    }
}

Vec2 TriangleMesh::get_surface_normal(index_t cell_index,
                                      index_t surface_id) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    Vec2 normal;
    switch (surface_id) {
        case 0:
            normal = Vec2(p2.y - p1.y, p1.x - p2.x);
            break;
        case 1:
            normal = Vec2(p3.y - p2.y, p2.x - p3.x);
            break;
        case 2:
            normal = Vec2(p1.y - p3.y, p3.x - p1.x);
            break;
        default:
            throw std::runtime_error("Invalid surface id");
    }

    return normal / normal.norm();
}

}  // namespace eigrene