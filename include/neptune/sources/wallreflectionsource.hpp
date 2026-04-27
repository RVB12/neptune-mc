#ifndef NEPTUNE_WALLREFLECTIONSOURCE_HPP
#define NEPTUNE_WALLREFLECTIONSOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/sources/source.hpp"

namespace neptune {

/**
 * @brief A source of particles that are generated from ions hitting a wall and
 * being reflected as neutrals.
 *
 * 1. The intensity of the source is determined by the flux of ions hitting the
 * wall, which is determined by integrating a drifting Maxwellian. A fraction Rt
 * of the ions will be reflected.
 *
 * 2. A fraction Rr of these reflected particles will be immediately
 * backscattered. To simulate these, we sample from the maxwellian of the ions,
 * but only consider particles that have a velocity towards the wall (this is a
 * sample from a truncated normal distribution). This velocity is then reflected
 * specularly to get the velocity of the backscattered particle. Along the way,
 * this velocity is also scaled to account for the electromagnetic sheath
 *
 * 3. The remaining fraction of the reflected particles will be released as
 * molecules. We assume that these immediately dissociate into neutrals with a
 * maxwellian distribution with a fixed temperature.
 *
 */
class WallReflectionSource : public Source {
   public:
    WallReflectionSource(index_t num_particles_,
                         index_t cell_idx_,
                         index_t surface_id_,
                         real_t sheath_potential_,
                         real_t franck_condon_temp_,
                         real_t Rt_,
                         real_t Rr_,
                         Species background_species,
                         Species particle_species,
                         const Mesh& mesh,
                         const Background& background);

    /**
     * @brief Generate a particle
     */
    Particle generate_particle(uint64_t particle_id,
                               rng_key_t key,
                               const Mesh& mesh,
                               const Background& background);

    index_t cell_idx = 0;    // The cell index where the particles are generated
    index_t surface_id = 0;  // The surface id where the particles are generated
    Species background_species;  // Species of the background particles that hit
                                 // the wall
    Species particle_species;    // Species of the generated particles
    real_t sheath_potential = 3.1;
    real_t franck_condon_temp =
        2.0;    // temperature of dissociated molecules in eV
    real_t Rt;  // fraction of ions that are reflected
    real_t Rr;  // fraction of reflected particles that are backscattered
};

}  // namespace neptune

#endif