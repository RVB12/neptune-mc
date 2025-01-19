#ifndef EIRON2_VOXELMESH_HPP
#define EIRON2_VOXELMESH_HPP

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vec2.hpp"
#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/geometry/particlegeometry.hpp"

namespace eiron2 {

class Index2D {
   public:
    Index2D(){};
    Index2D(index_t i_, index_t j_) : i(i_), j(j_) {}

    inline index_t& operator[](index_t ind)
    {
        assert(ind >= 0);
        assert(ind < 2);
        switch (ind) {
            case 0:
                return i;
            case 1:
                return j;
        }
    }

    // Member variables
    index_t i = -1;
    index_t j = -1;
};

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
    VoxelMesh(real_t Lx_, real_t Ly_, index_t nx_, index_t ny_)
        : Mesh(), Lx(Lx_), Ly(Ly_), nx(nx_), ny(ny_), dx(Lx / nx), dy(Ly / ny)
    {}

    index_t n_cells() const;

    CellCrossing get_next_crossing(ParticleGeometry& p) const;

    /**
     * @brief Convert a linear index to triplet index
     */
    Index2D linear_to_tuple(index_t index) const;

    /**
     * @brief Convert a triplet index to linear index
     */
    index_t tuple_to_linear(Index2D ij) const;

    /**
     * @brief Return the id of the cell that a given position is in.
     */
    index_t get_cell_index(Vec2 position) const;

    /**
     * @brief Return the size of the cell in each direction
     */
    Vec2 cell_size() const;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    bool is_in_cell(Vec2 p, index_t cell_index) const;

    /**
     * @brief Return the volume of a given cell
     */
    real_t get_volume(index_t cell_index) const;

    /**
     * @brief Return the normal of a given surface
     */
    Vec2 get_surface_normal(index_t cell_index, index_t surface_id) const;

    // Member variables
    real_t Lx = 0.;  // Length of the mesh in x direction
    real_t Ly = 0.;  // Length of the mesh in y direction
    index_t nx = 0;  // Number of cells in the x direction
    index_t ny = 0;  // Number of cells in the y direction
    real_t dx = 0.;  // Size of a cell in x direction
    real_t dy = 0.;  // Size of a cell in y direction
};

}  // namespace eiron2

#endif