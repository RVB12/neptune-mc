#ifndef NEPTUNE_POINTSOURCE_HPP
#define NEPTUNE_POINTSOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/sources/source.hpp"

namespace neptune {

/**
 * @brief A point source of particles
 *
 * The particles will be generated on a specified point.
 * The velocity is drawn from a normal distribution
 *
 */
class PointSource : public Source {
   public:
    PointSource(index_t num_particles_,
                real_t intensity_,
                Species species_,
                Vec3 pos_,
                Vec3 vel_mean_,
                real_t vel_stddev_)
        : Source(num_particles_),
          intensity(intensity_),
          species(species_),
          pos(pos_),
          vel_mean(vel_mean_),
          vel_stddev(vel_stddev_)
    {}

    /**
     * @brief Generate a particle
     */
    Particle generate_particle(uint64_t particle_id,
                               uint64_t seed,
                               const Mesh& mesh,
                               const Background& background);

    Species species;    // Species of the generated particles
    Vec3 pos;           // Initial position of the particles
    Vec3 vel_mean;      // Mean velocity of the particles
    real_t vel_stddev;  // Standard deviation of the velocity of the particles
    real_t intensity;   // The intensity of the source (particles per second)
    index_t cell_idx = -1;  // The cell index where the particles are generated
                            // initially set to -1 and updated to the true cell
                            // index during the particle generation process
};

}  // namespace neptune

#endif