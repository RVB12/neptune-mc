#include "eiron2/geometry/voxelmesh.hpp"

namespace eiron2 {

int VoxelMesh::n_cells() const { return nx * ny * nz; }

Index3D VoxelMesh::linear_to_triplet(int index) const
{
    int k = index / (nx * ny);
    int remainder = index % (nx * ny);
    int j = remainder / nx;
    int i = remainder % nx;

    return Index3D(i, j, k);
}

int VoxelMesh::triplet_to_linear(Index3D index) const
{
    return index.i + index.j * nx + index.k * nx * ny;
}

int VoxelMesh::find_cell(Vec3 p) const
{
    Index3D ind;
    assert(p.x >= 0 && p.x <= Lx);
    assert(p.y >= 0 && p.y <= Ly);
    assert(p.z >= 0 && p.z <= Lz);

    ind.i = p.x / dx;
    ind.j = p.y / dy;
    ind.k = p.z / dz;

    return triplet_to_linear(ind);
}

CellCrossing VoxelMesh::get_next_crossing(Particle& p) const
{
    Index3D ind = linear_to_triplet(p.cell_id);

    //
    // Calculate the time until collision for each surface
    //

    double tx;  // Time till collision with boundary in x direction
    if (p.velocity.x >= 0) {
        tx = ((ind.i + 1) * dx - p.position.x) / p.velocity.x;
    }
    else {
        tx = ((ind.i) * dx - p.position.x) / p.velocity.x;
    }

    double ty;  // Time till collision with boundary in y direction
    if (p.velocity.y >= 0) {
        ty = ((ind.j + 1) * dy - p.position.y) / p.velocity.y;
    }
    else {
        ty = ((ind.j) * dy - p.position.y) / p.velocity.y;
    }

    double tz;  // Time till collision with boundary in z direction
    if (p.velocity.z >= 0) {
        tz = ((ind.k + 1) * dy - p.position.y) / p.velocity.z;
    }
    else {
        tz = ((ind.k) * dy - p.position.y) / p.velocity.z;
    }

    //
    // Calculate the next cell index
    //

    Index3D next_ind = ind;
    double t;
    // Note that because of the < checks, it is not possible
    // for multiple of these branches to occur
    if (tx <= ty and tx <= tz) {
        t = tx;
        next_ind.i += (p.velocity.x >= 0) ? 1 : -1;
    }
    if (ty < tx and ty <= tz) {
        t = ty;
        next_ind.j += (p.velocity.y >= 0) ? 1 : -1;
    }
    if (tz < tx and tz < ty) {
        t = tz;
        next_ind.k += (p.velocity.z >= 0) ? 1 : -1;
    }

    // Calculate new position
    Vec3 pos_new = p.position + p.velocity * t;

    return CellCrossing(pos_new, triplet_to_linear(next_ind), t);

    // TODO: handle out of bounds
}

}  // namespace eiron2