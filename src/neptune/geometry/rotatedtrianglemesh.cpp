#include "neptune/geometry/rotatedtrianglemesh.hpp"

namespace neptune {

index_t RotatedTriangleMesh::n_cells() const { return triangles.size(); }

Vec2 carthesian_to_section(Vec3 position)
{
    real_t r = std::sqrt(position.x * position.x + position.y * position.y);
    return Vec2(r, position.z);
}

TriangleIndices RotatedTriangleMesh::get_triangle_indices(
    index_t cell_index) const
{
    return triangles[cell_index];
}

bool RotatedTriangleMesh::is_in_cell(Vec3 pos, index_t cell_index) const
{
    TriangleIndices tind = triangles[cell_index];
    Triangle t(vertices[tind.ind1], vertices[tind.ind2], vertices[tind.ind3]);

    return t.contains(carthesian_to_section(pos));
}

index_t RotatedTriangleMesh::get_cell_index(Vec3 pos) const
{
    Vec2 pos2d = carthesian_to_section(pos);

    auto candidates = spatialindex.query(pos2d);

    for (index_t i = 0; i < candidates.size(); i++) {
        index_t cell_index = candidates[i];
        if (is_in_cell(pos, cell_index)) {
            return cell_index;
        }
    }

    return -1;
}

real_t RotatedTriangleMesh::get_volume(index_t cell_index) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    // Area of the triangle in 2D
    // This calculation assumes the triangle is stored counter-clockwise
    real_t area_triangle = 0.5 * (-p2.y * p3.x + p1.y * (-p2.x + p3.x) +
                                  p1.x * (p2.y - p3.y) + p2.x * p3.y);

    // radial coordinate of the centroid of the triangle
    real_t centroid_r = (p1.x + p2.x + p3.x) / 3;

    // Volume of the rotated triangle
    // Pappus' centroid theorem
    return 2 * M_PI * centroid_r * area_triangle;
}

// Return the time until the intersection of particle with position pos and
// velocity vel with a cone defined by rotating points p1 and p2 around the
// z-axis. Note, this function does not verify whether the returned time is
// positive or whether the intersection lies between the two points. This should
// be done later.
real_t find_cone_intersection(const Vec2& p1,
                              const Vec2& p2,
                              const Vec3& pos,
                              const Vec3& vel,
                              bool is_on_surface)
{
    // Precompute slope terms
    real_t r_a = p1.x;
    real_t z_a = p1.y;
    real_t r_b = p2.x;
    real_t z_b = p2.y;
    real_t k_r = r_b - r_a;
    real_t k_z = z_b - z_a;

    if (std::abs(k_z) < 1.0e-10) {
        // The cone is actually a plane perpendicular to the z axis
        // Solve for t in the equation pos.z + t * vel.z = z_a
        if (is_on_surface or std::abs(vel.z) < 1.0e-12) {
            return -1;
        }
        real_t t = (z_a - pos.z) / vel.z;
        return t > 0 ? t : -1;
    }

    // Compute the coefficients of the quadratic equation
    real_t a =
        k_z * k_z * (vel.x * vel.x + vel.y * vel.y) - k_r * k_r * vel.z * vel.z;
    real_t b = 2 * k_z * k_z * (pos.x * vel.x + pos.y * vel.y) -
               2 * k_r * vel.z * (k_z * r_a + k_r * (pos.z - z_a));
    real_t c = k_z * k_z * (pos.x * pos.x + pos.y * pos.y) -
               std::pow(k_z * r_a + k_r * (pos.z - z_a), 2);

    // Solve the quadratic equation
    real_t discriminant = b * b - 4 * a * c;
    real_t t1;
    real_t t2;
    if (std::abs(discriminant) < 1.0e-16) {
        discriminant = 0.;
    }
    if (discriminant < 0.) {
        // No intersection, return negative number because it will be ignored.
        return -1;
    }
    else if (std::abs(a) < 1.0e-10) {
        // Only one solution
        // Use a different formula to avoid rounding errors
        t1 = -c / b;
        t2 = -c / b;
    }
    else if (is_on_surface) {
        // The initial position is on the cone (within tolerance)
        // One solution is at the initial position
        // Use a different formula to avoid rounding errors
        t1 = -b / a;
        t2 = 0;
    }
    else {
        // Normal case, just calculate both solutions
        if (b >= 0.) {
            t1 = (-b - std::sqrt(discriminant)) / (2 * a);
            t2 = c / (a * t1);
        }
        else {
            t1 = 2 * c / (-b + std::sqrt(discriminant));
            t2 = c / (a * t1);
        }
    }

    // Make sure t1 and t2 are not outside of the line segment
    real_t z1 = pos.z + t1 * vel.z;
    real_t zmin = std::min(z_a, z_b) - 1.0e-10;
    real_t zmax = std::max(z_a, z_b) + 1.0e-10;
    if (z1 < zmin or z1 > zmax) {
        t1 = -1;
    }
    real_t z2 = pos.z + t2 * vel.z;
    if (z2 < zmin or z2 > zmax) {
        t2 = -1;
    }

    // Return the smallest positive solution
    if (t1 > 0. and t2 > 0.) {
        return std::min(t1, t2);
    }
    else if (t1 > 0.) {
        return t1;
    }
    else if (t2 > 0.) {
        return t2;
    }
    else {
        return -1;
    }
}

CellCrossing RotatedTriangleMesh::get_next_crossing(ParticleGeometry& p) const
{
    assert(is_in_cell(p.position, p.cell_id));
// If the particle is on a boundary, make sure the velocity points inwards
#ifndef NDEBUG
    if (p.boundary_id == 0) {
        Vec3 normal = get_surface_normal(p.cell_id, 0, p.position);
        assert(normal.dot(p.velocity) <= 0);
    }
    else if (p.boundary_id == 1) {
        Vec3 normal = get_surface_normal(p.cell_id, 1, p.position);
        assert(normal.dot(p.velocity) <= 0);
    }
    else if (p.boundary_id == 2) {
        Vec3 normal = get_surface_normal(p.cell_id, 2, p.position);
        assert(normal.dot(p.velocity) <= 0);
    }
#endif

    TriangleIndices t = triangles[p.cell_id];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    real_t t0 = find_cone_intersection(p1, p2, p.position, p.velocity,
                                       p.boundary_id == 0);

    real_t t1 = find_cone_intersection(p2, p3, p.position, p.velocity,
                                       p.boundary_id == 1);

    real_t t2 = find_cone_intersection(p3, p1, p.position, p.velocity,
                                       p.boundary_id == 2);

    // t0, t1, and t2 are now the times until the intersection with the cones
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

    if (id_new < 0) {
        // If the id is negative, it is a boundary,
        // so negate it and subtract 1 (zero based) to get the boundary id
        assert(is_in_cell(pos_new, p.cell_id));
        crossing.boundary_condition_id = -id_new - 1;
        return crossing;
    }
    assert(is_in_cell(pos_new, id_new));

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

std::pair<index_t, index_t> RotatedTriangleMesh::get_boundary_crossing(
    Vec3 position, Vec3 direction) const
{
    // This function is not implemented yet
    throw std::runtime_error("get_boundary_crossing is not implemented");
}

Vec3 RotatedTriangleMesh::get_surface_normal(index_t cell_index,
                                             index_t surface_id,
                                             Vec3 position) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    real_t r = std::sqrt(position.x * position.x + position.y * position.y);
    real_t c = position.x / r;
    real_t s = position.y / r;

    Vec3 normal;
    real_t rr;
    switch (surface_id) {
        case 0:
            rr = p2.y - p1.y;
            normal = Vec3(rr * c, rr * s, p1.x - p2.x);
            break;
        case 1:
            rr = p3.y - p2.y;
            normal = Vec3(rr * c, rr * s, p2.x - p3.x);
            break;
        case 2:
            rr = p1.y - p3.y;
            normal = Vec3(rr * c, rr * s, p3.x - p1.x);
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

index_t RotatedTriangleMesh::get_n_surfaces(index_t cell_index) const
{
    // For triangles, there are always 3 surfaces
    return 3;
}

index_t RotatedTriangleMesh::get_neighbour(index_t cell_index,
                                           index_t surface_id) const
{
    assert(surface_id < 3);
    assert(cell_index < triangles.size());
    assert(surface_id >= 0);
    assert(cell_index >= 0);

    return neighbors[3 * cell_index + surface_id];
}

real_t RotatedTriangleMesh::get_surface_area(index_t cell_index,
                                             index_t surface_id) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    real_t r1, r2, z1, z2;
    switch (surface_id) {
        case 0:
            r1 = p1.x;
            r2 = p2.x;
            z1 = p1.y;
            z2 = p2.y;
            break;
        case 1:
            r1 = p2.x;
            r2 = p3.x;
            z1 = p2.y;
            z2 = p3.y;
            break;
        case 2:
            r1 = p3.x;
            r2 = p1.x;
            z1 = p3.y;
            z2 = p1.y;
            break;
        default:
            throw std::runtime_error("Invalid surface id");
    }

    // distance between the two points
    real_t s = std::sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));

    // Surface area of the cone frustum (https://mathworld.wolfram.com/ConicalFrustum.html)
    return M_PI * (r1 + r2) * s;
}

Vec3 RotatedTriangleMesh::get_center(index_t cell_index) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    return Vec3((p1.x + p2.x + p3.x) / 3, 0, (p1.y + p2.y + p3.y) / 3);
}

Vec3 RotatedTriangleMesh::get_surface_center(index_t cell_index,
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

Vec3 RotatedTriangleMesh::generate_random_position(index_t cell_index,
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

    real_t r = lambda1 * p1.x + lambda2 * p2.x + lambda3 * p3.x;
    real_t z = lambda1 * p1.y + lambda2 * p2.y + lambda3 * p3.y;

    // Generate random angle
    real_t theta = dist(rng) * 2 * M_PI;
    real_t x = r * std::cos(theta);
    real_t y = r * std::sin(theta);

    return Vec3(x, y, z);
}

Vec3 RotatedTriangleMesh::generate_random_surface_position(index_t cell_index,
                                                           index_t surface_id,
                                                           rng_t& rng) const
{
    TriangleIndices t = triangles[cell_index];
    Vec2 p1 = vertices[t.ind1];
    Vec2 p2 = vertices[t.ind2];
    Vec2 p3 = vertices[t.ind3];

    std::uniform_real_distribution<real_t> dist(0, 1);
    real_t u = dist(rng);

    // Generate a random angle along the surface
    real_t theta = dist(rng) * 2 * M_PI;
    real_t c = std::cos(theta);
    real_t s = std::sin(theta);

    if (surface_id == 0) {
        real_t r = p1.x + u * (p2.x - p1.x);
        return Vec3(c * r, s * r, p1.y + u * (p2.y - p1.y));
    }
    else if (surface_id == 1) {
        real_t r = p2.x + u * (p3.x - p2.x);
        return Vec3(c * r, s * r, p2.y + u * (p3.y - p2.y));
    }
    else if (surface_id == 2) {
        real_t r = p3.x + u * (p1.x - p3.x);
        return Vec3(c * r, s * r, p3.y + u * (p1.y - p3.y));
    }
    else {
        throw std::runtime_error("Invalid surface id");
    }
}

}  // namespace neptune