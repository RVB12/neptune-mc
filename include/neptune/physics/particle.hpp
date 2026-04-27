#ifndef NEPTUNE_PARTICLE_HPP
#define NEPTUNE_PARTICLE_HPP

#include <cassert>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/particlegeometry.hpp"
#include "neptune/geometry/vec3.hpp"
#include "neptune/physics/species.hpp"
#include "neptune/simulator/simulationparams.hpp"
namespace neptune {

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
     *            it is used to initialize the counter in the rng,
     *            if two particles have the same id, they will have the same rng
     * sequence, which will lead to incorrect results.
     */
    // initialize rng with high and low bits of the id for the first two parts
    // of the counter, and the key for the key
    Particle(uint64_t id_, rng_key_t key, Species species_)
        : id(id_),
          rng({{static_cast<rng_ctr_t::value_type>(id_ >> 32),
                static_cast<rng_ctr_t::value_type>(id_), 0, 0}},
              key),
          species(species_)
    {}

    // Member variables
    ParticleGeometry geometry;
    uint64_t id = 0;               // (unique) id of the particle
    real_t t = 0;                  // Time since start of simulation
    real_t weight = 1.0;           // Weight (used for continuous ionisation)
    Species species;  // Species of the particle
    real_t intensity = 1.0;  // How many particles this particle represents
    rng_t rng;               // Random number generator
    real_t u = 0;  // Scaled time till next background collision, this is
                   // typically drawn from an exp(1) distribution, but can be
                   // scaled for ARMC simulations
};

}  // namespace neptune

#endif