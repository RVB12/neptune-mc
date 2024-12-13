#ifndef EIRON2_MESH_HPP
#define EIRON2_MESH_HPP

#include "eiron2/geometry/cellcrossing.hpp"
#include "eiron2/geometry/vec3.hpp"
#include "eiron2/physics/particle.hpp"

namespace eiron2 {

/**
 * @brief This is a very general class for a mesh
 *
 */
class Mesh {
   public:
    Mesh() = default;
    ~Mesh() = default;

    /**
     * @brief Return the number of cells in the mesh
     */
    virtual int n_cells() const = 0;

    /**
     * For a given particle, if it continues at the current velocity,
     * find the position of the crossing with a cell boundary
     *
     * @param p
     */
    virtual CellCrossing get_next_crossing(Particle& p) const = 0;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    virtual bool is_in_cell(Vec3 p, int cell_index) const = 0;

    /**
     * Find the cell index corresponding to a given position.
     *
     * @param p
     */
    // virtual int get_cell_index(Vec3 position) const = 0;
};

}  // namespace eiron2

#endif