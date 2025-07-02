#include "neptune/geometry/voxelmesh.hpp"

namespace neptune {

index_t VoxelMesh::n_cells() const { return nx * ny * nz; }

Index3D VoxelMesh::linear_to_tuple(index_t index) const
{
    index_t k = index / (nx * ny);
    index_t rem = index % (nx * ny);
    index_t j = rem / nx;
    index_t i = rem % nx;

    return Index3D(i, j, k);
}

index_t VoxelMesh::tuple_to_linear(Index3D index) const
{
    // return index.i + index.j * nx;
    return index.i + index.j * nx + index.k * nx * ny;
}

bool VoxelMesh::is_in_cell(Vec3 p, index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);

    if (p.x < ind.i * dx || p.x > (ind.i + 1) * dx) return false;
    if (p.y < ind.j * dy || p.y > (ind.j + 1) * dy) return false;
    if (p.z < ind.k * dz || p.z > (ind.k + 1) * dz) return false;

    return true;
}

Vec3 VoxelMesh::cell_size() const { return Vec3(dx, dy, dz); }

real_t VoxelMesh::get_volume(index_t) const { return dx * dy * dz; }

index_t VoxelMesh::get_cell_index(Vec3 p) const
{
    Index3D ind;
    if (p.x >= Lx || p.y >= Ly || p.z >= Lz || p.x < 0 || p.y < 0 || p.z < 0) {
        return -1;
    }

    ind.i = p.x / dx;
    ind.j = p.y / dy;
    ind.k = p.z / dz;

    return tuple_to_linear(ind);
}

CellCrossing VoxelMesh::get_next_crossing(ParticleGeometry& p) const
{
    Index3D ind = linear_to_tuple(p.cell_id);

    //
    // Calculate the time until collision for each surface
    //

    real_t tx;  // Time till collision with boundary in x direction
    if (p.velocity.x >= 0) {
        tx = ((ind.i + 1) * dx - p.position.x) / p.velocity.x;
    }
    else {
        tx = ((ind.i) * dx - p.position.x) / p.velocity.x;
    }

    real_t ty;  // Time till collision with boundary in y direction
    if (p.velocity.y >= 0) {
        ty = ((ind.j + 1) * dy - p.position.y) / p.velocity.y;
    }
    else {
        ty = ((ind.j) * dy - p.position.y) / p.velocity.y;
    }

    real_t tz;  // Time till collision with boundary in z direction
    if (p.velocity.z >= 0) {
        tz = ((ind.k + 1) * dz - p.position.z) / p.velocity.z;
    }
    else {
        tz = ((ind.k) * dz - p.position.z) / p.velocity.z;
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
    if (next_ind.i >= nx) {
        crossing.boundary_condition_id = 3;
        return crossing;
    }
    if (next_ind.j >= ny) {
        crossing.boundary_condition_id = 4;
        return crossing;
    }
    if (next_ind.k >= nz) {
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

Vec3 VoxelMesh::get_surface_normal(index_t cell_index, index_t surface_id) const
{
    switch (surface_id) {
        case 0:
            return Vec3(1, 0, 0);
        case 1:
            return Vec3(0, 1, 0);
        case 2:
            return Vec3(0, 0, 1);
        case 3:
            return Vec3(-1, 0, 0);
        case 4:
            return Vec3(0, -1, 0);
        case 5:
            return Vec3(0, 0, -1);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

index_t VoxelMesh::get_n_surfaces(index_t) const
{
    // For voxels, there are always 6 surfaces
    return 6;
}


index_t VoxelMesh::get_neighbour(index_t cell_index, index_t surface_id) const{
    Index3D ind = linear_to_tuple(cell_index);

    switch (surface_id) {
        case 0:  // x == 0
            if (ind.i > 0) {
                return tuple_to_linear(Index3D(ind.i - 1, ind.j, ind.k));
            } else {
                return -1;
            }
        case 1:  // y == 0
            if (ind.j > 0) {
                return tuple_to_linear(Index3D(ind.i, ind.j - 1, ind.k));
            } else {
                return -2;
            }
        case 2:  // z == 0
            if (ind.k > 0) {
                return tuple_to_linear(Index3D(ind.i, ind.j, ind.k - 1));
            } else {
                return -3;
            }
        case 3:  // x == Lx
            if (ind.i < nx - 1) {
                return tuple_to_linear(Index3D(ind.i + 1, ind.j, ind.k));
            } else {
                return -4;
            }
        case 4:  // y == Ly
            if (ind.j < ny - 1) {
                return tuple_to_linear(Index3D(ind.i, ind.j + 1, ind.k));
            } else {
                return -5;
            }
        case 5:  // z == Lz
            if (ind.k < nz - 1) {
                return tuple_to_linear(Index3D(ind.i, ind.j, ind.k + 1));
            } else {
                return -6;
            }
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

real_t VoxelMesh::get_surface_area(index_t, index_t surface_id) const {

    switch (surface_id) {
        case 0:
        case 3:
            return dy * dz;
        case 1:
        case 4:
            return dx * dz;
        case 2:
        case 5:
            return dx * dy;
        default:
            throw std::runtime_error("Invalid surface id");
    }

}

Vec3 VoxelMesh::get_center(index_t cell_index) const
{
    Index3D ind = linear_to_tuple(cell_index);
    return Vec3((ind.i + 0.5) * dx, (ind.j + 0.5) * dy, (ind.k + 0.5) * dz);
}

Vec3 VoxelMesh::get_surface_center(index_t cell_index, index_t surface_id) const
{
    Index3D ind = linear_to_tuple(cell_index);
    switch (surface_id) {
        case 0:  // x == 0
            return Vec3(ind.i * dx, (ind.j + 0.5) * dy, (ind.k + 0.5) * dz);
        case 1:  // y == 0
            return Vec3((ind.i + 0.5) * dx, ind.j * dy, (ind.k + 0.5) * dz);
        case 2:  // z == 0
            return Vec3((ind.i + 0.5) * dx, (ind.j + 0.5) * dy, ind.k * dz);
        case 3:  // x == Lx
            return Vec3((ind.i + 1) * dx, (ind.j + 0.5) * dy, (ind.k + 0.5) * dz);
        case 4:  // y == Ly
            return Vec3((ind.i + 0.5) * dx, (ind.j + 1) * dy, (ind.k + 0.5) * dz);
        case 5:  // z == Lz
            return Vec3((ind.i + 0.5) * dx, (ind.j + 0.5) * dy, (ind.k + 1) * dz);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

std::pair<index_t, index_t> VoxelMesh::get_boundary_crossing(Vec3 position, Vec3 direction) const{
    //
    // Calculate the "time" until collision for each surface
    //

    real_t tx;  // Time till collision with boundary in x direction
    if (direction.x >= 0) {
        tx = (Lx - position.x) / direction.x;
    }
    else {
        tx = (-position.x) / direction.x;
    }

    real_t ty;  // Time till collision with boundary in y direction
    if (direction.y >= 0) {
        ty = (Ly - position.y) / direction.y;
    }
    else {
        ty = (-position.y) / direction.y;
    }

    real_t tz;  // Time till collision with boundary in z direction
    if (direction.z >= 0) {
        tz = (Lz - position.z) / direction.z;
    }
    else {
        tz = (-position.z) / direction.z;
    }

    // Calculate position of the boundary crossing
    Vec3 new_pos;
    index_t surface;
    if (tx <= ty and tx <= tz) {
        new_pos = position + direction * tx;
        surface = (direction.x >= 0) ? 3 : 0;
    }
    else if (ty <= tx and ty <= tz) {
        new_pos = position + direction * ty;
        surface = (direction.y >= 0) ? 4 : 1;
    }
    else {
        new_pos = position + direction * tz;
        surface = (direction.z >= 0) ? 5 : 2;
    }

    // Calculate index of the cell that crossing is in
    index_t cell_id = get_cell_index(new_pos);

    return std::make_pair(cell_id, surface);
}



}  // namespace neptune