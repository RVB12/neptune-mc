#ifndef NEPTUNE_VOXELMESH_HPP
#define NEPTUNE_VOXELMESH_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/cellcrossing.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/particlegeometry.hpp"
#include "neptune/geometry/index3D.hpp"

namespace neptune {

/**
 * @brief A structured rectangular mesh consisting of rectangular voxels
 *
 * The mesh is of size (Lx x Ly x Lz) and is split into (nx x ny x nz) cells.
 * The cells are ordered sequentially in the x, then the y and then the z
 * direction.
 *
 * The boundaries are numbered as:
 * 0 -> x == 0
 * 1 -> y == 0
 * 2 -> z == 0
 * 3 -> x == Lx
 * 4 -> y == Ly
 * 5 -> z == Lz
 *
 */
class VoxelMesh : public Mesh {
   public:
    VoxelMesh() = default;
    VoxelMesh(real_t Lx_, real_t Ly_, real_t Lz_, index_t nx_, index_t ny_, index_t nz_)
        : Mesh(), Lx(Lx_), Ly(Ly_), Lz(Lz_), nx(nx_), ny(ny_), nz(nz_), dx(Lx / nx), dy(Ly / ny), dz(Lz / nz)
    {}

    index_t n_cells() const;

    CellCrossing get_next_crossing(ParticleGeometry& p) const;

    /**
     * @brief Convert a linear index to triplet index
     */
    Index3D linear_to_tuple(index_t index) const;

    /**
     * @brief Convert a triplet index to linear index
     */
    index_t tuple_to_linear(Index3D ij) const;

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
    Vec3 get_surface_normal(index_t cell_index, index_t surface_id) const;

    /**
     * @brief Return the number of surfaces of a given cell
     */
    index_t get_n_surfaces(index_t cell_index) const;

    /**
     * @brief Return the index of the neighbour cell that shares a given surface
     */
    index_t get_neighbour(index_t cell_index, index_t surface_id) const;

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
     * @brief Return the first boundary crossing (i.e. an actual boundary, not a cell boundary)
     * That a particle will cross.
     */
    std::pair<index_t, index_t> get_boundary_crossing(Vec3 position, Vec3 direction) const;

    /**
     * @brief Generate a random position within a given cell
     */
    Vec3 generate_random_position(index_t cell_index, rng_t& rng) const;

    /**
     * @brief Generate a random position on a given surface of a cell
     */
    Vec3 generate_random_surface_position(index_t cell_index, index_t surface_id, rng_t& rng) const;

    // Member variables
    real_t Lx = 0.;  // Length of the mesh in x direction
    real_t Ly = 0.;  // Length of the mesh in y direction
    real_t Lz = 0.;  // Length of the mesh in z direction
    index_t nx = 0;  // Number of cells in the x direction
    index_t ny = 0;  // Number of cells in the y direction
    index_t nz = 0;  // Number of cells in the z direction
    real_t dx = 0.;  // Size of a cell in x direction
    real_t dy = 0.;  // Size of a cell in y direction
    real_t dz = 0.;  // Size of a cell in z direction
};

}  // namespace neptune

#endif