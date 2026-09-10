#include "neptune/geometry/rotatednonuniformmesh.hpp"

namespace neptune {

index_t RotatedNonUniformMesh::n_cells() const
{
    return (cb_r.size() - 1) * (cb_z.size() - 1);
}

Index3D RotatedNonUniformMesh::linear_to_tuple(index_t index) const
{
    index_t j = index / (cb_r.size() - 1);
    index_t i = index % (cb_r.size() - 1);

    return Index3D(i, j, 0);
}

index_t RotatedNonUniformMesh::tuple_to_linear(Index3D index) const
{
    return index.i + index.j * (cb_r.size() - 1);
}

bool RotatedNonUniformMesh::is_in_cell(Vec3 p, index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);

    real_t r = std::sqrt(p.x * p.x + p.y * p.y);
    real_t z = p.z;

    if (r < cb_r[ind.i] || r > cb_r[ind.i + 1]) return false;
    if (z < cb_z[ind.j] || z > cb_z[ind.j + 1]) return false;

    return true;
}

Vec3 RotatedNonUniformMesh::cell_size() const
{
    return Vec3(cb_r[1] - cb_r[0], cb_z[1] - cb_z[0], 0);
}

real_t RotatedNonUniformMesh::get_volume(index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);
    return M_PI * (cb_z[ind.k + 1] - cb_z[ind.k]) *
           (cb_r[ind.i + 1] * cb_r[ind.i + 1] - cb_r[ind.i] * cb_r[ind.i]);
}

index_t RotatedNonUniformMesh::get_cell_index(Vec3 p) const
{
    Index3D ind;

    real_t r = std::sqrt(p.x * p.x + p.y * p.y);
    real_t z = p.z;

    if (r >= cb_r[cb_r.size() - 1] || z >= cb_z[cb_z.size() - 1] || r < 0 ||
        z < cb_z[0]) {
        return -1;
    }

    // Use binary search to find the cell index
    ind.i = binary_search(cb_r, r);
    ind.j = binary_search(cb_z, z);
    ind.k = 0;

    return tuple_to_linear(ind);
}

/**
 * @brief Calculate the first positive intersection of a ray with a cylinder of
 * radius r_cyl centered on the z axis.
 *
 * @param p The position of the ray
 * @param v The direction of the ray
 * @param r_cyl The radius of the cylinder
 * @param is_on_surface Whether the ray is currently on the surface of the
 * cylinder
 *
 * @return real_t The time until the intersection, or -1 if there is no
 * intersection
 */
real_t ray_cylinder_intersection(Vec3 p,
                                 Vec3 v,
                                 real_t r_cyl,
                                 bool is_on_surface)
{
    real_t a = v.x * v.x + v.y * v.y;
    real_t b = 2 * (p.x * v.x + p.y * v.y);
    real_t c = p.x * p.x + p.y * p.y - r_cyl * r_cyl;

    real_t discriminant = b * b - 4 * a * c;

    real_t t1, t2;
    if (std::abs(discriminant) < 1e-14) {
        // If the discriminant is very small, we can treat it as zero
        t1 = t2 = -b / (2 * a);
    }
    else if (discriminant < 0) {
        // No real roots, no intersection
        return -1;
    }
    else if (std::abs(a) < 1e-14) {
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

CellCrossing RotatedNonUniformMesh::get_next_crossing(ParticleGeometry& p) const
{
    Index3D ind = linear_to_tuple(p.cell_id);

    //
    // Calculate the time until collision for each surface
    //

    // First, the easy part: the z surfaces
    real_t tz;  // Time till collision with boundary in z direction
    if (p.velocity.z >= 0) {
        tz = (cb_z[ind.k + 1] - p.position.z) / p.velocity.z;
    }
    else {
        tz = (cb_z[ind.k] - p.position.z) / p.velocity.z;
    }

    // Now, the r surfaces. Since these are curved, they form a hollow cylinder.
    // We need to solve the equation for the intersection of a line and a
    // cylinder to find the time until collision.
    real_t tr_inner = ray_cylinder_intersection(
        p.position, p.velocity, cb_r[ind.i], p.boundary_id == 0);
    real_t tr_outer = ray_cylinder_intersection(
        p.position, p.velocity, cb_r[ind.i + 1], p.boundary_id == 2);

    Index3D next_ind = ind;
    real_t t = 0;
    index_t surface;
    if (tr_inner > 0 and (tr_inner <= tr_outer or tr_outer < 0) and
        (tr_inner <= tz or tz < 0)) {
        t = tr_inner;
        next_ind.i -= 1;
        surface = 0;
    }
    else if (tr_outer > 0 and (tr_outer <= tr_inner or tr_inner < 0) and
             (tr_outer <= tz or tz < 0)) {
        t = tr_outer;
        next_ind.i += 1;
        surface = 2;
    }
    else if (tz > 0 and (tz <= tr_inner or tr_inner < 0) and
             (tz <= tr_outer or tr_outer < 0)) {
        t = tz;
        next_ind.j += (p.velocity.z >= 0) ? 1 : -1;
        surface = (p.velocity.z >= 0) ? 3 : 1;
    }
    else {
        // No valid collision found, this should not happen if the particle is
        // inside the cell
        throw std::runtime_error("No valid collision found");
    }

    CellCrossing crossing(p.position + p.velocity * t, t, surface);

    // Check for boundary crossings
    if (next_ind.i < 0) {
        crossing.boundary_condition_id = 0;
        return crossing;
    }
    if (next_ind.j < 0) {
        crossing.boundary_condition_id = 1;
        return crossing;
    }
    if (next_ind.i >= cb_r.size() - 1) {
        crossing.boundary_condition_id = 2;
        return crossing;
    }
    if (next_ind.j >= cb_z.size() - 1) {
        crossing.boundary_condition_id = 3;
        return crossing;
    }

    crossing.next_cell_id = tuple_to_linear(next_ind);
    return crossing;
}

Vec3 RotatedNonUniformMesh::get_surface_normal(index_t cell_index,
                                               index_t surface_id,
                                               Vec3 position) const
{
    Vec3 normal;
    switch (surface_id) {
        case 0:  // Inner r surface
            normal = Vec3(position.x, position.y, 0);
            normal = normal / normal.norm();
            break;
        case 1:  // lower z surface
            normal = Vec3(0, 0, 1);
            break;
        case 2:  // outer r surface
            normal = Vec3(-position.x, -position.y, 0);
            normal = normal / normal.norm();
            break;
        case 3:  // upper z surface
            normal = Vec3(0, 0, -1);
            break;
        default:
            throw std::runtime_error("Invalid surface id");
    }
    return normal;
}

index_t RotatedNonUniformMesh::get_n_surfaces(index_t) const
{
    // For the rectangles, there are always 4 surfaces
    return 4;
}

index_t RotatedNonUniformMesh::get_neighbour(index_t cell_index,
                                             index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);

    switch (surface_id) {
        case 0:  // r == inner
            if (ind.i > 0) {
                return tuple_to_linear(Index3D(ind.i - 1, ind.j, ind.k));
            }
            else {
                return -1;
            }
        case 1:  // z == lower
            if (ind.j > 0) {
                return tuple_to_linear(Index3D(ind.i, ind.j - 1, ind.k));
            }
            else {
                return -2;
            }
        case 2:  // r == outer
            if (ind.i < cb_r.size() - 2) {
                return tuple_to_linear(Index3D(ind.i + 1, ind.j, ind.k));
            }
            else {
                return -3;
            }
        case 3:  // z == upper
            if (ind.j < cb_z.size() - 2) {
                return tuple_to_linear(Index3D(ind.i, ind.j + 1, ind.k));
            }
            else {
                return -4;
            }
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

index_t RotatedNonUniformMesh::get_boundary_condition(index_t cell_index,
                                                     index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);

    switch (surface_id) {
        case 0:  // r == inner
            if (ind.i > 0) {
                return -1;
            }
            else {
                return 0;
            }
        case 1:  // z == lower
            if (ind.j > 0) {
                return -1;
            }
            else {
                return 1;
            }
        case 2:  // r == outer
            if (ind.i < cb_r.size() - 2) {
                return -1;
            }
            else {
                return 2;
            }
        case 3:  // z == upper
            if (ind.j < cb_z.size() - 2) {
                return -1;
            }
            else {
                return 3;
            }
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

real_t RotatedNonUniformMesh::get_surface_area(index_t cell_index,
                                               index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);

    // TODO: finish this code
    throw std::runtime_error("not implemented");
}

Vec3 RotatedNonUniformMesh::get_center(index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);
    return Vec3(0.5 * (cb_r[ind.i] + cb_r[ind.i + 1]), 0,
                0.5 * (cb_z[ind.k] + cb_z[ind.k + 1]));
}

Vec3 RotatedNonUniformMesh::get_surface_center(index_t cell_index,
                                               index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);

    // TODO: finish this code
    throw std::runtime_error("not implemented");
}

std::pair<index_t, index_t> RotatedNonUniformMesh::get_boundary_crossing(
    Vec3 position, Vec3 direction) const
{
    // This function is not implemented yet
    throw std::runtime_error("get_boundary_crossing is not implemented");
}

Vec3 RotatedNonUniformMesh::generate_random_position(index_t cell_index,
                                                     rng_t& rng) const
{
    Index3D ind = linear_to_tuple(cell_index);

    std::uniform_real_distribution<real_t> dist_r(cb_r[ind.i], cb_r[ind.i + 1]);
    std::uniform_real_distribution<real_t> dist_z(cb_z[ind.j], cb_z[ind.j + 1]);
    std::uniform_real_distribution<real_t> dist_phi(0, 2 * M_PI);

    real_t r = dist_r(rng);
    real_t phi = dist_phi(rng);
    real_t z = dist_z(rng);

    return Vec3(r * std::cos(phi), r * std::sin(phi), z);
}

Vec3 RotatedNonUniformMesh::generate_random_surface_position(index_t cell_index,
                                                             index_t surface_id,
                                                             rng_t& rng) const
{
    Index3D ind = linear_to_tuple(cell_index);

    std::uniform_real_distribution<real_t> dist_r(cb_r[ind.i], cb_r[ind.i + 1]);
    std::uniform_real_distribution<real_t> dist_z(cb_z[ind.j], cb_z[ind.j + 1]);
    std::uniform_real_distribution<real_t> dist_phi(0, 2 * M_PI);

    real_t r, phi, z;
    phi = dist_phi(rng);
    if (surface_id == 0) {
        r = cb_r[ind.i];
        z = dist_z(rng);
    }
    else if (surface_id == 1) {
        r = dist_r(rng);
        z = cb_z[ind.j];
    }
    else if (surface_id == 2) {
        r = cb_r[ind.i + 1];
        z = dist_z(rng);
    }
    else if (surface_id == 3) {
        r = dist_r(rng);
        z = cb_z[ind.j + 1];
    }
    else {
        throw std::runtime_error("Invalid surface id");
    }

    return Vec3(r * std::cos(phi), r * std::sin(phi), z);
}

}  // namespace neptune