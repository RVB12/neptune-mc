#ifndef NEPTUNE_SURFACESOURCE_HPP
#define NEPTUNE_SURFACESOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/sources/source.hpp"

namespace neptune {

/**
 * @brief A surface source of particles
 *
 * The particles will be generated somewhere on a specified surface of a
 * specified cell.
 *
 * The particle gets launched with a set velocity magnitude (v_norm) and a
 * direction sampled from a cosine distribution around the surface normal.
 *
 */
class SurfaceSource : public Source {
   public:
    SurfaceSource(index_t num_particles_,
                  real_t intensity_,
                  index_t cell_idx_,
                  index_t surface_id_,
                  real_t v_norm_)
        : Source(num_particles_),
          intensity(intensity_),
          cell_idx(cell_idx_),
          surface_id(surface_id_),
          v_norm(v_norm_)
    {}

    /**
     * @brief Generate a particle
     */
    Particle generate_particle(uint64_t particle_id,
                               rng_key_t key,
                               const Mesh& mesh,
                               const Background& background);

    real_t v_norm;         // Velocity magnitude of the particles
    real_t intensity;      // The intensity of the source (particles per second)
    index_t cell_idx = 0;  // The cell index where the particles are generated
    index_t surface_id = 0;  // The surface id where the particles are generated
};

}  // namespace neptune

#endif