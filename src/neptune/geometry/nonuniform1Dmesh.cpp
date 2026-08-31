#include "neptune/geometry/nonuniform1Dmesh.hpp"

namespace neptune {

index_t NonUniform1DMesh::n_cells() const { return (cb_x.size() - 1); }

Vec3 NonUniform1DMesh::cell_size() const
{
    return Vec3(cb_x[1] - cb_x[0], 1.0, 1.0);
}

bool NonUniform1DMesh::is_in_cell(Vec3 p, index_t cell_index) const
{
    if (p.x < cb_x[cell_index] || p.x > cb_x[cell_index + 1]) return false;
    return true;
}

real_t NonUniform1DMesh::get_volume(index_t cell_index) const
{
    return (cb_x[cell_index + 1] - cb_x[cell_index]);
}

index_t NonUniform1DMesh::get_cell_index(Vec3 p) const
{
    if (p.x >= cb_x[cb_x.size() - 1] || p.x < 0) {
        return -1;
    }

    return binary_search(cb_x, p.x);
}

CellCrossing NonUniform1DMesh::get_next_crossing(ParticleGeometry& p) const
{
    index_t ind = p.cell_id;

    //
    // Calculate the time until collision for each surface
    //

    real_t t;  // Time till collision with boundary in x direction
    index_t surface;
    index_t next_ind;
    if (p.velocity.x >= 0) {
        t = (cb_x[ind + 1] - p.position.x) / p.velocity.x;
        surface = 1;
        next_ind = ind + 1;
    }
    else {
        t = (cb_x[ind] - p.position.x) / p.velocity.x;
        surface = 0;
        next_ind = ind - 1;
    }

    // Calculate new position
    Vec3 pos_new = p.position + p.velocity * t;

    CellCrossing crossing(pos_new, t, surface);

    if (next_ind < 0) {
        crossing.next_cell_id = -1;
        crossing.boundary_condition_id = 0;
        return crossing;
    }
    if (next_ind >= cb_x.size() - 1) {
        crossing.next_cell_id = -1;
        crossing.boundary_condition_id = 1;
        return crossing;
    }

    crossing.next_cell_id = next_ind;

    // If the current surface crossing is on the left, that means
    // it is the right surface of the next cell, and vice versa
    if (surface == 0) {
        crossing.next_surface_id = 1;
    }
    else {
        crossing.next_surface_id = 0;
    }

    return crossing;
}

Vec3 NonUniform1DMesh::get_surface_normal(index_t cell_index,
                                          index_t surface_id,
                                          Vec3 position) const
{
    switch (surface_id) {
        case 0:  // x == 0
            return Vec3(1, 0, 0);
        case 1:  // x == Lx
            return Vec3(-1, 0, 0);
        default:
            throw std::runtime_error("Invalid surface id");
    }
}

index_t NonUniform1DMesh::get_n_surfaces(index_t) const { return 2; }

index_t NonUniform1DMesh::get_neighbour(index_t cell_index,
                                        index_t surface_id) const
{
    if (surface_id == 0) {
        if (cell_index == 0) {
            return -1;
        }
        else {
            return cell_index - 1;
        }
    }
    if (surface_id == 1) {
        if (cell_index == cb_x.size() - 2) {
            return -1;
        }
        else {
            return cell_index + 1;
        }
    }
    throw std::runtime_error("Invalid surface id");
}

real_t NonUniform1DMesh::get_surface_area(index_t cell_index,
                                          index_t surface_id) const
{
    return 1.0;  // The surface area is always 1 in the y and z directions
}

Vec3 NonUniform1DMesh::get_center(index_t cell_index) const
{
    return Vec3(0.5 * (cb_x[cell_index] + cb_x[cell_index + 1]), 0.0,
                0.0);  // The center in y and z is always 0
}

Vec3 NonUniform1DMesh::get_surface_center(index_t cell_index,
                                          index_t surface_id) const
{
    return Vec3(cb_x[cell_index + surface_id], 0.0, 0.0);
}

std::pair<index_t, index_t> NonUniform1DMesh::get_boundary_crossing(
    Vec3 position, Vec3 direction) const
{
    // This function is not implemented yet
    throw std::runtime_error("get_boundary_crossing is not implemented");
}

Vec3 NonUniform1DMesh::generate_random_position(index_t cell_index,
                                                rng_t& rng) const
{
    std::uniform_real_distribution<real_t> dist(0, 1);
    real_t x = cb_x[cell_index] + (dist(rng)) * (cb_x[cell_index + 1] - cb_x[cell_index]);
    return Vec3(x, 0.0, 0.0);
}

Vec3 NonUniform1DMesh::generate_random_surface_position(index_t cell_index,
                                                        index_t surface_id,
                                                        rng_t& rng) const
{
    return Vec3(cb_x[cell_index + surface_id], 0.0, 0.0);
}

}  // namespace neptune