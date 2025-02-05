#ifndef EIGRENE_MESH_HPP
#define EIGRENE_MESH_HPP

#include "eigrene/geometry/cellcrossing.hpp"
#include "eigrene/geometry/particlegeometry.hpp"
#include "eigrene/geometry/vec2.hpp"
#include "eigrene/common/utils.hpp"

namespace eigrene {

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
    virtual index_t n_cells() const = 0;

    /**
     * For a given particle, if it continues at the current velocity,
     * find the position of the crossing with a cell boundary
     *
     * @param p
     */
    virtual CellCrossing get_next_crossing(ParticleGeometry& p) const = 0;

    /**
     * Check if a given position is inside of the cell with the given index
     */
    virtual bool is_in_cell(Vec2 p, index_t cell_index) const = 0;

    /**
     * Find the cell index corresponding to a given position.
     * Note: this may be very slow for unstructured meshes.
     * 
     * returns -1 if no cell was found.
     */
    virtual index_t get_cell_index(Vec2 position) const = 0;

    /**
     * @brief Return the volume of a given cell
     */
    virtual real_t get_volume(index_t cell_index) const = 0;

    /**
     * @brief Return the normal of a given surface
     */
    virtual Vec2 get_surface_normal(index_t cell_index, index_t surface_id) const = 0;
};

}  // namespace eigrene

#endif