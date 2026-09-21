#include "neptune/geometry/trianglemesh.hpp"

namespace neptune {

index_t TriangleMesh::n_cells() const { return triangles.size(); }

TriangleIndices TriangleMesh::get_triangle_indices(index_t cell_index) const
{
    return triangles[cell_index];
}

bool TriangleMesh::is_in_cell(Vec3 pos, index_t cell_index) const
{
    TriangleIndices tind = triangles[cell_index];
    Triangle t(vertices[tind.ind1], vertices[tind.ind2], vertices[tind.ind3]);

    return t.contains(Vec2(pos.x, pos.z));
}

index_t TriangleMesh::get_cell_index(Vec3 pos) const
{
    Vec2 pos2d(pos.x, pos.z);

    auto candidates = spatialindex.query(pos2d);

    for (index_t i = 0; i < candidates.size(); i++) {
        index_t cell_index = candidates[i];
        if (is_in_cell(pos, cell_index)) {
            return cell_index;
        }
    }

    return -1;
}

void TriangleMesh::find_particle(ParticleGeometry& p) const
{
    index_t cell_index = get_cell_index(p.position);
    if (cell_index == -1) {
        throw std::runtime_error(
            "Particle is outside of the mesh. Cannot find cell.");
    }

    // Get distance to each edge of the triangle
    TriangleIndices tind = triangles[cell_index];
    Vec2 p1 = vertices[tind.ind1];
    Vec2 p2 = vertices[tind.ind2];
    Vec2 p3 = vertices[tind.ind3];

    Vec2 pos2d(p.position.x, p.position.z);

    real_t d0 = std::abs((p2.y - p1.y) * pos2d.x - (p2.x - p1.x) * pos2d.y +
                         p2.x * p1.y - p2.y * p1.x) /
                std::sqrt((p2.y - p1.y) * (p2.y - p1.y) +
                          (p2.x - p1.x) * (p2.x - p1.x));
    real_t d1 = std::abs((p3.y - p2.y) * pos2d.x - (p3.x - p2.x) * pos2d.y +
                         p3.x * p2.y - p3.y * p2.x) /
                std::sqrt((p3.y - p2.y) * (p3.y - p2.y) +
                          (p3.x - p2.x) * (p3.x - p2.x));
    real_t d2 = std::abs((p1.y - p3.y) * pos2d.x - (p1.x - p3.x) * pos2d.y +
                         p1.x * p3.y - p1.y * p3.x) /
                std::sqrt((p1.y - p3.y) * (p1.y - p3.y) +
                          (p1.x - p3.x) * (p1.x - p3.x));

    index_t boundary_id = -1;
    if (d0 < BOUNDARY_TOLERANCE) {
        boundary_id = 0;
    }
    else if (d1 < BOUNDARY_TOLERANCE) {
        boundary_id = 1;
    }
    else if (d2 < BOUNDARY_TOLERANCE) {
        boundary_id = 2;
    }

    if (boundary_id != -1) {
        // Check if the particle should be in the current cell or the neighbor
        // cell
        Vec3 normal = get_surface_normal(cell_index, boundary_id, p.position);
        // normal points into the cell, so if component in direction of normal
        // is negative, the particle is moving out of the cell
        if (normal.dot(p.velocity) < 0) {
            p.cell_id = neighbors[3 * cell_index + boundary_id];
            // Find boundary id corresponding to this edge
            for (index_t i = 0; i < 3; ++i) {
                if (neighbors[3 * p.cell_id + i] == cell_index) {
                    p.boundary_id = i;
                    break;
                }
            }
        }
        else {
            p.cell_id = cell_index;
            p.boundary_id = 0;
        }
    }
    else {
        p.boundary_id = -1;
    }
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
    // assert(is_in_cell(p.position, p.cell_id));
    TriangleIndices t = triangles[p.cell_id];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    Vec2 pos2d(p.position.x, p.position.z);
    Vec2 vel2d(p.velocity.x, p.velocity.z);

    // Time till intersection with edge 0
    real_t t0;
    if (p.boundary_id == 0) {
        t0 = -1;
    }
    else {
        t0 = find_edge_intersection(p1, p2, pos2d, vel2d);
    }

    real_t t1;
    if (p.boundary_id == 1) {
        t1 = -1;
    }
    else {
        t1 = find_edge_intersection(p2, p3, pos2d, vel2d);
    }

    real_t t2;
    if (p.boundary_id == 2) {
        t2 = -1;
    }
    else {
        t2 = find_edge_intersection(p3, p1, pos2d, vel2d);
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

    Vec3 pos_new = p.position + p.velocity * t_new;

    CellCrossing crossing = CellCrossing(pos_new, t_new, surface_id);

    crossing.boundary_condition_id =
        get_boundary_condition(p.cell_id, surface_id);

    if (id_new >= 0) {
        crossing.next_cell_id = id_new;
        // Find the surface id of the next cell
        for (index_t i = 0; i < 3; i++) {
            if (neighbors[3 * id_new + i] == p.cell_id) {
                crossing.next_surface_id = i;
                break;
            }
        }
    }

    return crossing;
}

std::pair<index_t, index_t> TriangleMesh::get_boundary_crossing(
    Vec3 position, Vec3 direction) const
{
    // This function is not implemented yet
    throw std::runtime_error("get_boundary_crossing is not implemented");
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

        for (index_t j = 0; j < candidates.size(); j++) {
            index_t k = candidates[j];
            if (i == k) continue;
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

Vec3 TriangleMesh::get_surface_normal(index_t cell_index,
                                      index_t surface_id,
                                      Vec3 position) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    Vec3 normal;
    switch (surface_id) {
        case 0:
            normal = Vec3(p2.y - p1.y, 0, p1.x - p2.x);
            break;
        case 1:
            normal = Vec3(p3.y - p2.y, 0, p2.x - p3.x);
            break;
        case 2:
            normal = Vec3(p1.y - p3.y, 0, p3.x - p1.x);
            break;
        default:
            throw std::runtime_error("Invalid surface id");
    }

    // Calculate the signed area
    // If this is negative, the calculated normal will point inwards
    // If this is positive, the calculated normal will point outwards, in that
    // case, we need to flip the normal
    real_t signed_area =
        (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);

    if (signed_area > 0) {
        normal = -normal;
    }

    return normal / normal.norm();
}

index_t TriangleMesh::get_n_surfaces(index_t cell_index) const
{
    // For triangles, there are always 3 surfaces
    return 3;
}

index_t TriangleMesh::get_neighbour(index_t cell_index,
                                    index_t surface_id) const
{
    assert(surface_id < 3);
    assert(cell_index < triangles.size());
    assert(surface_id >= 0);
    assert(cell_index >= 0);

    return neighbors[3 * cell_index + surface_id];
}

index_t TriangleMesh::get_boundary_condition(index_t cell_index,
                                             index_t surface_id) const
{
    assert(surface_id < 3);
    assert(cell_index < triangles.size());
    assert(surface_id >= 0);
    assert(cell_index >= 0);

    return boundary_conditions[3 * cell_index + surface_id];
}

real_t TriangleMesh::get_surface_area(index_t cell_index,
                                      index_t surface_id) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    if (surface_id == 0) {
        return (p2 - p1).norm();
    }
    else if (surface_id == 1) {
        return (p3 - p2).norm();
    }
    else if (surface_id == 2) {
        return (p1 - p3).norm();
    }
    else {
        throw std::runtime_error("Invalid surface id");
    }
}

Vec3 TriangleMesh::get_center(index_t cell_index) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    return Vec3((p1.x + p2.x + p3.x) / 3, 0, (p1.y + p2.y + p3.y) / 3);
}

Vec3 TriangleMesh::get_surface_center(index_t cell_index,
                                      index_t surface_id) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    switch (surface_id) {
        case 0:
            return Vec3((p1.x + p2.x) / 2, 0, (p1.y + p2.y) / 2);
        case 1:
            return Vec3((p2.x + p3.x) / 2, 0, (p2.y + p3.y) / 2);
        case 2:
            return Vec3((p3.x + p1.x) / 2, 0, (p3.y + p1.y) / 2);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

Vec3 TriangleMesh::generate_random_position(index_t cell_index,
                                            rng_t& rng) const
{
    std::uniform_real_distribution<real_t> dist(0, 1);

    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    // Generate random barycentric coordinates
    // This should result in nice uniform sampling within the triangle
    // (need this check this though)
    real_t u = dist(rng);
    real_t v = dist(rng);

    real_t s = std::sqrt(u);
    real_t lambda1 = 1.0 - s;
    real_t lambda2 = v * s;
    real_t lambda3 = 1.0 - lambda1 - lambda2;

    real_t x = lambda1 * p1.x + lambda2 * p2.x + lambda3 * p3.x;
    real_t z = lambda1 * p1.y + lambda2 * p2.y + lambda3 * p3.y;
    real_t y = 0;  // Since this is a 2D triangle mesh, we set y to 0

    return Vec3(x, y, z);
}

Vec3 TriangleMesh::generate_random_surface_position(index_t cell_index,
                                                    index_t surface_id,
                                                    rng_t& rng) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    std::uniform_real_distribution<real_t> dist(0, 1);
    real_t u = dist(rng);

    if (surface_id == 0) {
        return Vec3(p1.x + u * (p2.x - p1.x), 0, p1.y + u * (p2.y - p1.y));
    }
    else if (surface_id == 1) {
        return Vec3(p2.x + u * (p3.x - p2.x), 0, p2.y + u * (p3.y - p2.y));
    }
    else if (surface_id == 2) {
        return Vec3(p3.x + u * (p1.x - p3.x), 0, p3.y + u * (p1.y - p3.y));
    }
    else {
        throw std::runtime_error("Invalid surface id");
    }
}

}  // namespace neptune