#ifndef NEPTUNE_SURFACESOURCE_HPP
#define NEPTUNE_SURFACESOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/source.hpp"

namespace neptune {

/**
 * @brief A surface source of particles
 *
 * The particles will be generated somewhere on a specified surface of a
 * specified cell.
 *
 */
class SurfaceSource : public Source {
   public:
    SurfaceSource(index_t num_particles_,
                  real_t intensity_,
                  index_t cell_idx_,
                  index_t surface_id_,
                  Vec3 vel_mean_,
                  real_t vel_stddev_)
        : Source(num_particles_, intensity_),
          cell_idx(cell_idx_),
          surface_id(surface_id_),
          vel_mean(vel_mean_),
          vel_stddev(vel_stddev_)
    {}

    /**
     * @brief Generate a particle
     */
    Particle generate_particle(index_t i_particle, const Mesh& mesh);

    Vec3 vel_mean;      // Mean velocity of the particles
    real_t vel_stddev;  // Standard deviation of the velocity of the particles
    index_t cell_idx = 0;  // The cell index where the particles are generated
    index_t surface_id = 0;  // The surface id where the particles are generated
};

}  // namespace neptune

#endif