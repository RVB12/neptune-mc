#include "eiron2/geometry/trianglemesh.hpp"

namespace eiron2 {

int TriangleMesh::n_cells() const { return triangles.size(); }

Triangle TriangleMesh::get_triangle(int cell_index) const
{
    return triangles[cell_index];
}

inline double triangle_sign_helper(const Vec2& p1,
                                   const Vec2& p2,
                                   const Vec2& p3)
{
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

// Source:
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
bool TriangleMesh::is_in_cell(Vec2 pos, int cell_index) const
{
    Triangle t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    // Half-plane checks
    double d1 = triangle_sign_helper(pos, p1, p2);
    double d2 = triangle_sign_helper(pos, p2, p3);
    double d3 = triangle_sign_helper(pos, p3, p1);

    // Check if the point is on the same side of all triangle edges
    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

int TriangleMesh::get_cell_index(Vec2 pos) const
{
    // For now, this is really stupid, we're just manually going through all the
    // triangles and checking which one the position is in. A kd tree would
    // massively speed this up
    for (int i = 0; i < n_cells(); i++) {
        if (is_in_cell(pos, i)) {
            return i;
        }
    }

    return -1;
}

double TriangleMesh::get_volume(int cell_index) const
{
    Triangle t = triangles[cell_index];
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
inline double find_edge_intersection(const Vec2& p1,
                                     const Vec2& p2,
                                     const Vec2& pos,
                                     const Vec2& vel)
{
    // Construct 2x2 matrix equation to find the intersection
    double a00 = p2.x - p1.x;
    double a10 = p2.y - p1.y;
    double a01 = -vel.x;
    double a11 = -vel.y;

    double b0 = pos.x - p1.x;
    double b1 = pos.y - p1.y;

    double determinant = a00 * a11 - a01 * a10;

    if (std::abs(determinant) < 1.0e-8) {
        // matrix is singular, no intersection
        // return a negative number because those will be ignored anyway.
        return -1;
    }

    // solve system using cramer's rule.
    return (a00 * b1 - b0 * a10) / determinant;
}

CellCrossing TriangleMesh::get_next_crossing(ParticleGeometry& p) const
{
    Triangle t = triangles[p.cell_id];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    double t1 = find_edge_intersection(p1, p2, p.position, p.velocity);
    double t2 = find_edge_intersection(p2, p3, p.position, p.velocity);
    double t3 = find_edge_intersection(p3, p1, p.position, p.velocity);

    std::cout << "position is in cell? " << is_in_cell(p.position, p.cell_id) << std::endl;
    std::cout << "position: " << p.position << std::endl;
    std::cout << "velocity: " << p.velocity << std::endl;
    std::cout << "cell_id: " << p.cell_id << std::endl;
    std::cout << "found cell id: " << get_cell_index(p.position) << std::endl;
    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;
    std::cout << "t3: " << t3 << std::endl;

    double t_new;
    int id_new;
    if (t1 > 0 and t1 >= t2 and t1 >= t3) {
        t_new = t1;
        id_new = neighbors[3 * p.cell_id];
    }
    else if (t2 > 0 and t2 >= t1 and t2 >= t3) {
        t_new = t2;
        id_new = neighbors[3 * p.cell_id + 1];
    }
    else if (t3 > 0 and t3 >= t1 and t3 >= t2) {
        t_new = t3;
        id_new = neighbors[3 * p.cell_id + 2];
    }
    else {
        throw std::runtime_error("Failed to find next crossing");
    }

    Vec2 pos_new = p.position + p.velocity * t_new;

    int boundary_id = -1;
    if (id_new < 0) {
        boundary_id = -id_new;
        id_new = -1;
    }

    std::cout << "id_new: " << id_new << std::endl;

    return CellCrossing(pos_new, t_new, id_new, boundary_id);
}

void TriangleMesh::make_neighborlist()
{
    // For now, this is really stupid, we're just manually going through all the
    // triangles and checking the edges
    for (int i = 0; i < n_cells(); i++) {
        Triangle t1 = triangles[i];

        // Initialize to -1, both to keep track if all neighbors of this triangle have been found,
        // but also to keep track of boundaries.
        neighbors[3 * i] = -1;
        neighbors[3 * i + 1] = -1;
        neighbors[3 * i + 2] = -1;

        for (int j = 0; j < n_cells(); j++) {
            if (i == j) continue;

            Triangle t2 = triangles[j];

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

}  // namespace eiron2