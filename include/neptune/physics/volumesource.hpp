#ifndef NEPTUNE_VOLUMESOURCE_HPP
#define NEPTUNE_VOLUMESOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/source.hpp"

namespace neptune {

/**
 * @brief A volume source of particles
 *
 * The particles will be generated somewhere in a specified cell.
 */
class VolumeSource : public Source {
   public:
    VolumeSource(index_t num_particles_,
                 real_t intensity_,
                 index_t cell_idx_,
                 Vec3 vel_mean_,
                 real_t vel_stddev_)
        : Source(num_particles_),
          intensity(intensity_),
          cell_idx(cell_idx_),
          vel_mean(vel_mean_),
          vel_stddev(vel_stddev_)
    {}

    /**
     * @brief Generate a particle
     */
    Particle generate_particle(index_t i_particle, const Mesh& mesh, const Background& background);

    Vec3 vel_mean;      // Mean velocity of the particles
    real_t intensity;  // The intensity of the source (particles per second)
    real_t vel_stddev;  // Standard deviation of the velocity of the particles
    index_t cell_idx = 0;  // The cell index where the particles are generated
};

}  // namespace neptune

#endif