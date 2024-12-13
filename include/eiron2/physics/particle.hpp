#ifndef EIRON2_PARTICLE_HPP
#define EIRON2_PARTICLE_HPP

#include <cassert>

#include "eiron2/geometry/vec3.hpp"

namespace eiron2 {

/**
 * @brief Class containing basic info of a particle, position and velocity
 */
class Particle {
   public:
    // Contructors
    Particle() = default;
    Particle(Vec3 position_, Vec3 velocity_, int cell_id_)
        : position(position_), velocity(velocity_), cell_id(cell_id_)
    {}

    // Member variables
    Vec3 position;     // Current position
    Vec3 velocity;     // Current velocity
    int cell_id = -1;  // id of the cell the particle is in
    int id = -1;       // id of the particle
};

}  // namespace eiron2

#endif