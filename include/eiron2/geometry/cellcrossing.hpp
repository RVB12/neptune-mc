#ifndef EIRON2_CROSSING_HPP
#define EIRON2_CROSSING_HPP

#include <cassert>

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
    CellCrossing(Vec2 position_, double t_, int next_cell_id_, int boundary_id_)
        : position(position_),
          next_cell_id(next_cell_id_),
          boundary_id(boundary_id_),
          t(t_)
    {}

    // Member variables
    Vec2 position;  // The position of the crossing
    int next_cell_id =
        -1;  // the id of the next cell, -1 if it is a boundary crossing
    int boundary_id =
        -1;    // the id of the boundary, -1 if it is a cell crossing
    double t;  // Time until the crossing
};

}  // namespace eiron2

#endif