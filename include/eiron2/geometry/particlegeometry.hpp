#ifndef EIRON2_PARTICLEGEOMETRY_HPP
#define EIRON2_PARTICLEGEOMETRY_HPP

#include <cassert>

#include "eiron2/common/vec2.hpp"

namespace eiron2 {

/**
 * @brief Class containing geometric information about a particle,
 * this consists of the particle's position, its velocity and
 * the id of the cell the particle is in.
 */
class ParticleGeometry {
   public:
    // Contructors
    ParticleGeometry() = default;
    ParticleGeometry(Vec2 position_, Vec2 velocity_, int cell_id_)
        : position(position_), velocity(velocity_), cell_id(cell_id_)
    {}

    // Member variables
    Vec2 position;     // Current position
    Vec2 velocity;     // Current velocity
    int cell_id = -1;  // id of the cell the particle is in
};

}  // namespace eiron2

#endif