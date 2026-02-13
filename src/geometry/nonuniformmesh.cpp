#include "neptune/geometry/nonuniformmesh.hpp"

namespace neptune {

index_t NonUniformMesh::n_cells() const
{
    return (cb_x.size() - 1) * (cb_y.size() - 1) * (cb_z.size() - 1);
}

Index3D NonUniformMesh::linear_to_tuple(index_t index) const
{
    index_t k = index / (cb_x.size() - 1) / (cb_y.size() - 1);
    index_t rem = index % ((cb_x.size() - 1) * (cb_y.size() - 1));
    index_t j = rem / (cb_x.size() - 1);
    index_t i = rem % (cb_x.size() - 1);

    return Index3D(i, j, k);
}

index_t NonUniformMesh::tuple_to_linear(Index3D index) const
{
    // return index.i + index.j * nx;
    return index.i + index.j * (cb_y.size() - 1) +
           index.k * (cb_y.size() - 1) * (cb_z.size() - 1);
}

bool NonUniformMesh::is_in_cell(Vec3 p, index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);

    if (p.x < cb_x[ind.i] || p.x > cb_x[ind.i + 1]) return false;
    if (p.y < cb_y[ind.j] || p.y > cb_y[ind.j + 1]) return false;
    if (p.z < cb_z[ind.k] || p.z > cb_z[ind.k + 1]) return false;

    return true;
}

Vec3 NonUniformMesh::cell_size() const
{
    return Vec3(cb_x[1] - cb_x[0], cb_y[1] - cb_y[0], cb_z[1] - cb_z[0]);
}

real_t NonUniformMesh::get_volume(index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);
    return (cb_x[ind.i + 1] - cb_x[ind.i]) * (cb_y[ind.j + 1] - cb_y[ind.j]) *
           (cb_z[ind.k + 1] - cb_z[ind.k]);
}

index_t NonUniformMesh::get_cell_index(Vec3 p) const
{
    Index3D ind;
    if (p.x >= cb_x[cb_x.size() - 1] || p.y >= cb_y[cb_y.size() - 1] ||
        p.z >= cb_z[cb_z.size() - 1] || p.x < 0 || p.y < 0 || p.z < 0) {
        return -1;
    }

    // Use binary search to find the cell index
    ind.i = binary_search(cb_x, p.x);
    ind.j = binary_search(cb_y, p.y);
    ind.k = binary_search(cb_z, p.z);

    return tuple_to_linear(ind);
}

CellCrossing NonUniformMesh::get_next_crossing(ParticleGeometry& p) const
{
    Index3D ind = linear_to_tuple(p.cell_id);

    //
    // Calculate the time until collision for each surface
    //

    real_t tx;  // Time till collision with boundary in x direction
    if (p.velocity.x >= 0) {
        tx = (cb_x[ind.i + 1] - p.position.x) / p.velocity.x;
    }
    else {
        tx = (cb_x[ind.i] - p.position.x) / p.velocity.x;
    }

    real_t ty;  // Time till collision with boundary in y direction
    if (p.velocity.y >= 0) {
        ty = (cb_y[ind.j + 1] - p.position.y) / p.velocity.y;
    }
    else {
        ty = (cb_y[ind.j] - p.position.y) / p.velocity.y;
    }

    real_t tz;  // Time till collision with boundary in z direction
    if (p.velocity.z >= 0) {
        tz = (cb_z[ind.k + 1] - p.position.z) / p.velocity.z;
    }
    else {
        tz = (cb_z[ind.k] - p.position.z) / p.velocity.z;
    }

    //
    // Calculate the next cell index
    //

    Index3D next_ind = ind;
    real_t t = 0;
    index_t surface;
    if (tx <= ty and tx <= tz) {
        t = tx;
        next_ind.i += (p.velocity.x >= 0) ? 1 : -1;
        surface = (p.velocity.x >= 0) ? 3 : 0;
    }
    else if (ty <= tx and ty <= tz) {
        t = ty;
        next_ind.j += (p.velocity.y >= 0) ? 1 : -1;
        surface = (p.velocity.y >= 0) ? 4 : 1;
    }
    else {
        t = tz;
        next_ind.k += (p.velocity.z >= 0) ? 1 : -1;
        surface = (p.velocity.z >= 0) ? 5 : 2;
    }

    // Calculate new position
    Vec3 pos_new = p.position + p.velocity * t;

    CellCrossing crossing(pos_new, t, surface);

    // Check for boundary crossings
    if (next_ind.i < 0) {
        crossing.boundary_condition_id = 0;
        return crossing;
    }
    if (next_ind.j < 0) {
        crossing.boundary_condition_id = 1;
        return crossing;
    }
    if (next_ind.k < 0) {
        crossing.boundary_condition_id = 2;
        return crossing;
    }
    if (next_ind.i >= cb_x.size() - 1) {
        crossing.boundary_condition_id = 3;
        return crossing;
    }
    if (next_ind.j >= cb_y.size() - 1) {
        crossing.boundary_condition_id = 4;
        return crossing;
    }
    if (next_ind.k >= cb_z.size() - 1) {
        crossing.boundary_condition_id = 5;
        return crossing;
    }

    crossing.next_cell_id = tuple_to_linear(next_ind);
    switch (surface) {
        case 0:
            crossing.next_surface_id = 3;
            break;
        case 1:
            crossing.next_surface_id = 4;
            break;
        case 2:
            crossing.next_surface_id = 5;
            break;
        case 3:
            crossing.next_surface_id = 0;
            break;
        case 4:
            crossing.next_surface_id = 1;
            break;
        case 5:
            crossing.next_surface_id = 2;
            break;
    }

    return crossing;
}

Vec3 NonUniformMesh::get_surface_normal(index_t cell_index,
                                        index_t surface_id,
                                        Vec3 position) const
{
    switch (surface_id) {
        case 0: // x == 0
            return Vec3(1, 0, 0);
        case 1: // y == 0
            return Vec3(0, 1, 0);
        case 2: // z == 0
            return Vec3(0, 0, 1);
        case 3: // x == Lx
            return Vec3(-1, 0, 0);
        case 4: // y == Ly
            return Vec3(0, -1, 0);
        case 5: // z == Lz
            return Vec3(0, 0, -1);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

index_t NonUniformMesh::get_n_surfaces(index_t) const
{
    // For voxels, there are always 6 surfaces
    return 6;
}

index_t NonUniformMesh::get_neighbour(index_t cell_index,
                                      index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);

    switch (surface_id) {
        case 0:  // x == 0
            if (ind.i > 0) {
                return tuple_to_linear(Index3D(ind.i - 1, ind.j, ind.k));
            }
            else {
                return -1;
            }
        case 1:  // y == 0
            if (ind.j > 0) {
                return tuple_to_linear(Index3D(ind.i, ind.j - 1, ind.k));
            }
            else {
                return -2;
            }
        case 2:  // z == 0
            if (ind.k > 0) {
                return tuple_to_linear(Index3D(ind.i, ind.j, ind.k - 1));
            }
            else {
                return -3;
            }
        case 3:  // x == Lx
            if (ind.i < cb_x.size() - 2) {
                return tuple_to_linear(Index3D(ind.i + 1, ind.j, ind.k));
            }
            else {
                return -4;
            }
        case 4:  // y == Ly
            if (ind.j < cb_y.size() - 2) {
                return tuple_to_linear(Index3D(ind.i, ind.j + 1, ind.k));
            }
            else {
                return -5;
            }
        case 5:  // z == Lz
            if (ind.k < cb_z.size() - 2) {
                return tuple_to_linear(Index3D(ind.i, ind.j, ind.k + 1));
            }
            else {
                return -6;
            }
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

real_t NonUniformMesh::get_surface_area(index_t cell_index,
                                        index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);
    switch (surface_id) {
        case 0:
        case 3:
            return (cb_y[ind.j + 1] - cb_y[ind.j]) *
                   (cb_z[ind.k + 1] - cb_z[ind.k]);
        case 1:
        case 4:
            return (cb_x[ind.i + 1] - cb_x[ind.i]) *
                   (cb_z[ind.k + 1] - cb_z[ind.k]);
        case 2:
        case 5:
            return (cb_x[ind.i + 1] - cb_x[ind.i]) *
                   (cb_y[ind.j + 1] - cb_y[ind.j]);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

Vec3 NonUniformMesh::get_center(index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);
    return Vec3(0.5 * (cb_x[ind.i] + cb_x[ind.i + 1]),
                0.5 * (cb_y[ind.j] + cb_y[ind.j + 1]),
                0.5 * (cb_z[ind.k] + cb_z[ind.k + 1]));
}

Vec3 NonUniformMesh::get_surface_center(index_t cell_index,
                                        index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);
    switch (surface_id) {
        case 0:  // x == 0
            return Vec3(cb_x[ind.i], 0.5 * (cb_y[ind.j] + cb_y[ind.j + 1]),
                        0.5 * (cb_z[ind.k] + cb_z[ind.k + 1]));
        case 1:  // y == 0
            return Vec3(0.5 * (cb_x[ind.i] + cb_x[ind.i + 1]), cb_y[ind.j],
                        0.5 * (cb_z[ind.k] + cb_z[ind.k + 1]));
        case 2:  // z == 0
            return Vec3(0.5 * (cb_x[ind.i] + cb_x[ind.i + 1]),
                        0.5 * (cb_y[ind.j] + cb_y[ind.j + 1]), cb_z[ind.k]);
        case 3:  // x == Lx
            return Vec3(cb_x[ind.i + 1], 0.5 * (cb_y[ind.j] + cb_y[ind.j + 1]),
                        0.5 * (cb_z[ind.k] + cb_z[ind.k + 1]));
        case 4:  // y == Ly
            return Vec3(0.5 * (cb_x[ind.i] + cb_x[ind.i + 1]), cb_y[ind.j + 1],
                        0.5 * (cb_z[ind.k] + cb_z[ind.k + 1]));
        case 5:  // z == Lz
            return Vec3(0.5 * (cb_x[ind.i] + cb_x[ind.i + 1]),
                        0.5 * (cb_y[ind.j] + cb_y[ind.j + 1]), cb_z[ind.k + 1]);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

std::pair<index_t, index_t> NonUniformMesh::get_boundary_crossing(
    Vec3 position, Vec3 direction) const
{
    // This function is not implemented yet
    throw std::runtime_error("get_boundary_crossing is not implemented");
}

Vec3 NonUniformMesh::generate_random_position(index_t cell_index,
                                              rng_t& rng) const
{
    Index3D ind = linear_to_tuple(cell_index);

    std::uniform_real_distribution<real_t> dist(0, 1);

    real_t x = cb_x[ind.i] + (dist(rng)) * (cb_x[ind.i + 1] - cb_x[ind.i]);
    real_t y = cb_y[ind.j] + (dist(rng)) * (cb_y[ind.j + 1] - cb_y[ind.j]);
    real_t z = cb_z[ind.k] + (dist(rng)) * (cb_z[ind.k + 1] - cb_z[ind.k]);

    return Vec3(x, y, z);
}

Vec3 NonUniformMesh::generate_random_surface_position(index_t cell_index,
                                                      index_t surface_id,
                                                      rng_t& rng) const
{
    Index3D ind = linear_to_tuple(cell_index);

    std::uniform_real_distribution<real_t> dist(0, 1);

    switch (surface_id) {
        case 0:  // x == 0
            return Vec3(
                cb_x[ind.i],
                cb_y[ind.j] + (dist(rng)) * (cb_y[ind.j + 1] - cb_y[ind.j]),
                cb_z[ind.k] + (dist(rng)) * (cb_z[ind.k + 1] - cb_z[ind.k]));
        case 1:  // y == 0
            return Vec3(
                cb_x[ind.i] + (dist(rng)) * (cb_x[ind.i + 1] - cb_x[ind.i]),
                cb_y[ind.j],
                cb_z[ind.k] + (dist(rng)) * (cb_z[ind.k + 1] - cb_z[ind.k]));
        case 2:  // z == 0
            return Vec3(
                cb_x[ind.i] + (dist(rng)) * (cb_x[ind.i + 1] - cb_x[ind.i]),
                cb_y[ind.j] + (dist(rng)) * (cb_y[ind.j + 1] - cb_y[ind.j]),
                cb_z[ind.k]);
        case 3:  // x == Lx
            return Vec3(
                cb_x[ind.i + 1],
                cb_y[ind.j] + (dist(rng)) * (cb_y[ind.j + 1] - cb_y[ind.j]),
                cb_z[ind.k] + (dist(rng)) * (cb_z[ind.k + 1] - cb_z[ind.k]));
        case 4:  // y == Ly
            return Vec3(
                cb_x[ind.i] + (dist(rng)) * (cb_x[ind.i + 1] - cb_x[ind.i]),
                cb_y[ind.j + 1],
                cb_z[ind.k] + (dist(rng)) * (cb_z[ind.k + 1] - cb_z[ind.k]));
        case 5:  // z == Lz
            return Vec3(
                cb_x[ind.i] + (dist(rng)) * (cb_x[ind.i + 1] - cb_x[ind.i]),
                cb_y[ind.j] + (dist(rng)) * (cb_y[ind.j + 1] - cb_y[ind.j]),
                cb_z[ind.k + 1]);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

}  // namespace neptune