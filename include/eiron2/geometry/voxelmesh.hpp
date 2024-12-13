#ifndef EIRON2_VOXELMESH_HPP
#define EIRON2_VOXELMESH_HPP

#include <tuple>

#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/mesh.hpp"
#include "eiron2/geometry/vec3.hpp"
#include "eiron2/physics/particle.hpp"

namespace eiron2 {

class Index3D {
   public:
    Index3D(){};
    Index3D(int i_, int j_, int k_) : i(i_), j(j_), k(k_) {}

    inline int& operator[](int ind)
    {
        assert(ind >= 0);
        assert(ind < 3);
        switch (ind) {
            case 0:
                return i;
            case 1:
                return j;
            case 2:
                return k;
        }
    }

    // Member variables
    int i = -1;
    int j = -1;
    int k = -1;
};

/**
 * @brief A structured mesh consisting of rectangular voxels
 *
 */
class VoxelMesh : public Mesh {
   public:
    VoxelMesh(double Lx_, double Ly_, double Lz_, int nx_, int ny_, int nz_)
        : Mesh(),
          Lx(Lx_),
          Ly(Ly_),
          Lz(Lz_),
          nx(nx_),
          ny(ny_),
          nz(nz_),
          dx(Lx / nx),
          dy(Ly / ny),
          dz(Lz / nz)
    {}

    int n_cells() const;

    CellCrossing get_next_crossing(Particle& p) const;

    /**
     * @brief Convert a linear index to triplet index
     */
    Index3D linear_to_triplet(int index) const;

    /**
     * @brief Convert a triplet index to linear index
     */
    int triplet_to_linear(Index3D ijk) const;

    /**
     * @brief Return the id of the cell that a given position is in.
     */
    int find_cell(Vec3 position) const;

   private:
    // Member variables
    double Lx = 0.;  // Length of the mesh in x direction
    double Ly = 0.;  // Length of the mesh in y direction
    double Lz = 0.;  // Length of the mesh in z direction
    int nx = 0;      // Number of cells in the x direction
    int ny = 0;      // Number of cells in the y direction
    int nz = 0;      // Number of cells in the z direction
    double dx = 0.;  // Size of a cell in x direction
    double dy = 0.;  // Size of a cell in y direction
    double dz = 0.;  // Size of a cell in z direction
};

}  // namespace eiron2

#endif