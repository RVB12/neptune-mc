#ifndef NEPTUNE_NONUNIFORM1DMESH_HPP
#define NEPTUNE_NONUNIFORM1DMESH_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/cellcrossing.hpp"
#include "neptune/geometry/index3D.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/particlegeometry.hpp"
#include "neptune/geometry/vec2.hpp"

namespace neptune {

/**
 * @brief A structured rectangular mesh consisting of rectangular voxels
 *
 * The mesh is of size (Lx x 1 x 1) and is split into nx cells.
 * 
 * Each cell is only boundary on the left and the right side, and the mesh
 * is infinite in the y and z directions.
 * However, for the purpose of volume/density, the cell has length 1 in the y and z directions.
 * 
 * i.e. This is a 1D mesh and the y and z coordinates are ignored.
 *
 */
class NonUniform1DMesh : public Mesh {
   public:
    NonUniform1DMesh() = default;
    NonUniform1DMesh(vector<real_t>& cb_x)
        : Mesh(), cb_x(cb_x)
    {}

    NonUniform1DMesh(vector<real_t>&& cb_x)
        : Mesh(),
          cb_x(std::move(cb_x))
    {}

    NonUniform1DMesh(NonUniform1DMesh&& other) noexcept
        : Mesh(std::move(other)),
          cb_x(std::move(other.cb_x))
    {}

    NonUniform1DMesh(const NonUniform1DMesh& other) = default;

    NonUniform1DMesh& operator=(NonUniform1DMesh&& other) noexcept
    {
        if (this != &other) {
            Mesh::operator=(std::move(other));
            cb_x = std::move(other.cb_x);
        }
        return *this;
    }

    index_t n_cells() const;

    CellCrossing get_next_crossing(ParticleGeometry& p) const;

    /**
     * @brief Return the id of the cell that a given position is in.
     */
    index_t get_cell_index(Vec3 position) const;

    /**
     * @brief Return the size of the cell in each direction
     */
    Vec3 cell_size() const;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    bool is_in_cell(Vec3 p, index_t cell_index) const;

    /**
     * @brief Return the volume of a given cell
     */
    real_t get_volume(index_t cell_index) const;

    /**
     * @brief Return the normal of a given surface
     */
    Vec3 get_surface_normal(index_t cell_index,
                            index_t surface_id,
                            Vec3 position) const;

    /**
     * @brief Return the number of surfaces of a given cell
     */
    index_t get_n_surfaces(index_t cell_index) const;

    /**
     * @brief Return the index of the neighbour cell that shares a given surface
     */
    index_t get_neighbour(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the boundary condition of a given surface
     */
    index_t get_boundary_condition(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the area of a given surface
     */
    real_t get_surface_area(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the center of a given cell
     */
    Vec3 get_center(index_t cell_index) const;

    /**
     * @brief Return the center of a given surface of a cell
     */
    Vec3 get_surface_center(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the first boundary crossing (i.e. an actual boundary, not a
     * cell boundary) That a particle will cross.
     */
    std::pair<index_t, index_t> get_boundary_crossing(Vec3 position,
                                                      Vec3 direction) const;

    /**
     * @brief Generate a random position within a given cell
     */
    Vec3 generate_random_position(index_t cell_index, rng_t& rng) const;

    /**
     * @brief Generate a random position on a given surface of a cell
     */
    Vec3 generate_random_surface_position(index_t cell_index,
                                          index_t surface_id,
                                          rng_t& rng) const;

    /**
     * @brief Return whether the mesh uses toroidal coordinates or carthesian
     * coordinates
     */
    bool is_toroidal() const { return false; }

    // Member variables
    // Note, these vectors must be sorted!
    vector<real_t> cb_x;  // x-coordinates of the cell boundaries
};

}  // namespace neptune

#endif