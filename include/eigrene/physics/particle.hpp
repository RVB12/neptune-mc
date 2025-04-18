#ifndef EIGRENE_PARTICLE_HPP
#define EIGRENE_PARTICLE_HPP

#include <cassert>

#include "eigrene/common/utils.hpp"
#include "eigrene/geometry/vec3.hpp"
#include "eigrene/geometry/particlegeometry.hpp"

namespace eigrene {

/**
 * @brief Class containing all the information of a particle
 */
class Particle {
   public:
    /**
     * @brief Construct a new Particle object
     * 
     * @param id_ A unique id for the particle
     *            Note: it is super important that this id is unique.
     *            it is used to seed the random number generator for
     *            this particle so if two particles have the same id
     *            it will lead to hard to debug errors!!
     */
    Particle(index_t id_) : id(id_), rng(id) {}

    // Member variables
    ParticleGeometry geometry;
    index_t id = -1;      // (unique) id of the particle
    real_t t = 0;         // Time
    real_t weight = 1.0;  // Weight (used for continuous ionisation)
    real_t mass = 1.0;    // Mass of the particle [kg]
    real_t intensity = 1.0; // How many particles this particle represents
    rng_t rng;            // Random number generator
};

}  // namespace eigrene

#endif