#ifndef EIGRENE_CROSSING_HPP
#define EIGRENE_CROSSING_HPP

#include <cassert>

#include "eigrene/common/utils.hpp"
#include "eigrene/geometry/vec2.hpp"

namespace eigrene {

/**
 * @brief Class containing information about a cell crossing
 *        A cell crossing occurs when a particle moves from
 *        one cell to another. Or when the particles collides
 *        with the boundary of the grid.
 */
class CellCrossing {
   public:
    // Contructors
    CellCrossing() = default;
    CellCrossing(Vec2 position_, real_t t_, index_t surface_id_)
        : position(position_), t(t_), surface_id(surface_id_)
    {}
    CellCrossing(Vec2 position_,
                 real_t t_,
                 index_t surface_id_,
                 index_t next_cell_id_,
                 index_t next_surface_id_,
                 index_t boundary_condition_id_)
        : position(position_),
          t(t_),
          surface_id(surface_id_),
          next_cell_id(next_cell_id_),
          next_surface_id(next_surface_id_),
          boundary_condition_id(boundary_condition_id_)
    {}

    //
    // Member variables
    //

    // The position of the crossing
    Vec2 position;
    // The time it will take to reach the crossing
    real_t t = 0;
    // The id of the surface in the current cell the particle is crossing
    index_t surface_id = -1;
    // The id of the next cell
    // -1 if there is no next cell
    index_t next_cell_id = -1;
    // The id of the boundary the particle will be on in the next cell
    // -1 if there is no next cell
    index_t next_surface_id = -1;
    // The id of the boundary condition
    // -1 if this is a normal crossing between cells
    index_t boundary_condition_id = -1;
};

}  // namespace eigrene

#endif