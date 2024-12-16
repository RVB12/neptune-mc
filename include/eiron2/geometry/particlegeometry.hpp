#ifndef EIRON2_PARTICLEGEOMETRY_HPP
#define EIRON2_PARTICLEGEOMETRY_HPP

#include <cassert>

#include "eiron2/common/vec3.hpp"

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
    ParticleGeometry(Vec3 position_, Vec3 velocity_, int cell_id_)
        : position(position_), velocity(velocity_), cell_id(cell_id_)
    {}

    // Member variables
    Vec3 position;     // Current position
    Vec3 velocity;     // Current velocity
    int cell_id = -1;  // id of the cell the particle is in
};

}  // namespace eiron2

#endif