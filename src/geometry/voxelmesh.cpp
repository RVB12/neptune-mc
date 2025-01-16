#include "eiron2/geometry/voxelmesh.hpp"

namespace eiron2 {

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
    assert(p.x >= 0 && p.x <= Lx);
    assert(p.y >= 0 && p.y <= Ly);

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
    if (tx <= ty) {
        t = tx;
        next_ind.i += (p.velocity.x >= 0) ? 1 : -1;
    } else {
        t = ty;
        next_ind.j += (p.velocity.y >= 0) ? 1 : -1;
    }

    // Calculate new position
    Vec2 pos_new = p.position + p.velocity * t;

    // Check for boundary crossings
    if(next_ind.i < 0){
        return CellCrossing(pos_new, t, -1, 0);
    }
    if(next_ind.j < 0){
        return CellCrossing(pos_new, t, -1, 1);
    }
    if(next_ind.i >= nx){
        return CellCrossing(pos_new, t, -1, 2);
    }
    if(next_ind.j >= ny){
        return CellCrossing(pos_new, t, -1, 3);
    }

    return CellCrossing(pos_new, t, tuple_to_linear(next_ind), -1);
}

}  // namespace eiron2