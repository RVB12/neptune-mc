#ifndef EIRON2_CROSSING_HPP
#define EIRON2_CROSSING_HPP

#include <cassert>

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vec2.hpp"

namespace eiron2 {

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
    CellCrossing(Vec2 position_,
                 real_t t_,
                 index_t surface_id,
                 index_t next_cell_id_,
                 index_t boundary_id_)
        : position(position_),
          surface_id(surface_id),
          next_cell_id(next_cell_id_),
          boundary_id(boundary_id_),
          t(t_)
    {}

    // Member variables
    Vec2 position;  // The position of the crossing
    // index_t old_cell_id =
    //     -1;  // the id of the cell the particle is coming from
    index_t surface_id = -1;  // the id of the surface in the current cell the
                              // particle is crossing
    index_t next_cell_id =
        -1;  // the id of the next cell, -1 if it is a boundary crossing
    index_t boundary_id =
        -1;    // the id of the boundary, -1 if it is a cell crossing
    real_t t;  // Time until the crossing
};

}  // namespace eiron2

#endif