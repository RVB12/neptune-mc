#include "eigrene/geometry/voxelmesh.hpp"

namespace eigrene {

index_t VoxelMesh::n_cells() const { return nx * ny; }

Index2D VoxelMesh::linear_to_tuple(index_t index) const
{
    // index_t j = index / nx;
    // index_t i = index % nx;

    index_t i = index / ny;
    index_t j = index % ny;

    return Index2D(i, j);
}

index_t VoxelMesh::tuple_to_linear(Index2D index) const
{
    // return index.i + index.j * nx;
    return index.i * ny + index.j;
}

bool VoxelMesh::is_in_cell(Vec2 p, index_t cell_index) const{
    Index2D ind = linear_to_tuple(cell_index);

    if(p.x < ind.i * dx || p.x > (ind.i + 1) * dx) return false;
    if(p.y < ind.j * dy || p.y > (ind.j + 1) * dy) return false;

    return true;
}

Vec2 VoxelMesh::cell_size() const { return Vec2(dx, dy); }

real_t VoxelMesh::get_volume(index_t) const {
    return dx * dy;
}

index_t VoxelMesh::get_cell_index(Vec2 p) const
{
    Index2D ind;
    if(p.x >= Lx || p.y >= Ly || p.x < 0 || p.y < 0){
        return -1;
    }

    ind.i = p.x / dx;
    ind.j = p.y / dy;

    return tuple_to_linear(ind);
}

CellCrossing VoxelMesh::get_next_crossing(ParticleGeometry& p) const
{
    Index2D ind = linear_to_tuple(p.cell_id);

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

    //
    // Calculate the next cell index
    //

    Index2D next_ind = ind;
    real_t t = 0;
    index_t surface;
    if (tx <= ty) {
        t = tx;
        next_ind.i += (p.velocity.x >= 0) ? 1 : -1;
        surface = (p.velocity.x >= 0) ? 2 : 0;
    } else {
        t = ty;
        next_ind.j += (p.velocity.y >= 0) ? 1 : -1;
        surface = (p.velocity.y >= 0) ? 3 : 1;
    }

    // Calculate new position
    Vec2 pos_new = p.position + p.velocity * t;

    CellCrossing crossing(pos_new, t, surface);

    // Check for boundary crossings
    if(next_ind.i < 0){
        crossing.boundary_condition_id = 0;
        return crossing;
    }
    if(next_ind.j < 0){
        crossing.boundary_condition_id = 1;
        return crossing;
    }
    if(next_ind.i >= nx){
        crossing.boundary_condition_id = 2;
        return crossing;
    }
    if(next_ind.j >= ny){
        crossing.boundary_condition_id = 3;
        return crossing;
    }

    crossing.next_cell_id = tuple_to_linear(next_ind);
    if(surface == 0){
        crossing.next_surface_id = 2;
    } else if (surface == 1){
        crossing.next_surface_id = 3;
    } else if (surface == 2){
        crossing.next_surface_id = 0;
    } else if (surface == 3){
        crossing.next_surface_id = 1;
    }

    return crossing;
}

Vec2 VoxelMesh::get_surface_normal(index_t cell_index, index_t surface_id) const
{
    switch (surface_id) {
        case 0:
            return Vec2(-1, 0);
        case 1:
            return Vec2(0, -1);
        case 2:
            return Vec2(1, 0);
        case 3:
            return Vec2(0, 1);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

}  // namespace eigrene