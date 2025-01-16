#ifndef EIRON2_PARTICLE_HPP
#define EIRON2_PARTICLE_HPP

#include <cassert>

#include "eiron2/common/utils.hpp"
#include "eiron2/geometry/particlegeometry.hpp"
#include "eiron2/common/vec3.hpp"

namespace eiron2 {

/**
 * @brief Class containing all the information of a particle
 */
class Particle {
   public:
    // Contructors
    Particle() = default;
    Particle(ParticleGeometry geometry_) : geometry(geometry_) {}

    // Member variables
    ParticleGeometry geometry;
    index_t id = -1;          // (unique) id of the particle
    real_t t = 0;         // Time
    real_t weight = 1.0;  // Weight (used for continuous ionisation)
    real_t mass = 1.0;    // Mass of the particle [kg]
};

}  // namespace eiron2

#endif