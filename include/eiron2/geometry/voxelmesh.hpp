#ifndef EIRON2_VOXELMESH_HPP
#define EIRON2_VOXELMESH_HPP

#include <tuple>

#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/geometry/particlegeometry.hpp"
#include "eiron2/common/vec2.hpp"

namespace eiron2 {

class Index2D {
   public:
    Index2D(){};
    Index2D(int i_, int j_) : i(i_), j(j_) {}

    inline int& operator[](int ind)
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
    int i = -1;
    int j = -1;
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
    VoxelMesh(double Lx_, double Ly_, int nx_, int ny_)
        : Mesh(),
          Lx(Lx_),
          Ly(Ly_),
          nx(nx_),
          ny(ny_),
          dx(Lx / nx),
          dy(Ly / ny)
    {}

    int n_cells() const;

    CellCrossing get_next_crossing(ParticleGeometry& p) const;

    /**
     * @brief Convert a linear index to triplet index
     */
    Index2D linear_to_tuple(int index) const;

    /**
     * @brief Convert a triplet index to linear index
     */
    int tuple_to_linear(Index2D ij) const;

    /**
     * @brief Return the id of the cell that a given position is in.
     */
    int get_cell_index(Vec2 position) const;

    /**
     * @brief Return the size of the cell in each direction
     */
    Vec2 cell_size() const;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    bool is_in_cell(Vec2 p, int cell_index) const;

    /**
     * @brief Return the volume of a given cell
     */
    double get_volume(int cell_index) const;

   private:
    // Member variables
    double Lx = 0.;  // Length of the mesh in x direction
    double Ly = 0.;  // Length of the mesh in y direction
    int nx = 0;      // Number of cells in the x direction
    int ny = 0;      // Number of cells in the y direction
    double dx = 0.;  // Size of a cell in x direction
    double dy = 0.;  // Size of a cell in y direction
};

}  // namespace eiron2

#endif