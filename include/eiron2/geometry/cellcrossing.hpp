#ifndef EIRON2_CROSSING_HPP
#define EIRON2_CROSSING_HPP

#include <cassert>

#include "vec3.hpp"

namespace eiron2 {

/**
 * @brief Class containing information about a cell crossing
 *        A cell crossing occurs when a particle moves from
 *        one cell to another.
 */
class CellCrossing {
   public:
    // Contructors
    CellCrossing() = default;
    CellCrossing(Vec3 position_, int next_cell_id_, double t_)
        : position(position_), next_cell_id(next_cell_id_), t(t_)
    {}

    // Member variables
    Vec3 position;          // The position of the crossing
    int next_cell_id = -1;  // the id of the next cell
    double t;               // Time until the crossing
};

}  // namespace eiron2

#endif