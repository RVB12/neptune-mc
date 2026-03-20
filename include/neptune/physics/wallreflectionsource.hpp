#ifndef NEPTUNE_WALLREFLECTIONSOURCE_HPP
#define NEPTUNE_WALLREFLECTIONSOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/source.hpp"
#include "neptune/physics/physicalconstants.hpp"

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
                         real_t ion_vel_x_,
                         real_t ion_vel_y_,
                         real_t ion_vel_z_,
                         real_t ion_temp_,
                         real_t ion_density_,
                         real_t sheath_potential_,
                         real_t franck_condon_temp_,
                         real_t Rt_,
                         real_t Rr_)
        : Source(num_particles_, 0),
          cell_idx(cell_idx_),
          surface_id(surface_id_),
          ion_vel_x(ion_vel_x_),
          ion_vel_y(ion_vel_y_),
          ion_vel_z(ion_vel_z_),
          ion_temp(ion_temp_),
          ion_density(ion_density_),
          sheath_potential(sheath_potential_),
          franck_condon_temp(franck_condon_temp_),
          Rt(Rt_),
          Rr(Rr_)
    {}

    /**
     * @brief Generate a particle
     */
    Particle generate_particle(index_t i_particle, const Mesh& mesh);

    index_t cell_idx = 0;    // The cell index where the particles are generated
    index_t surface_id = 0;  // The surface id where the particles are generated
    real_t ion_vel_x;        // ion velocity in m/s
    real_t ion_vel_y;
    real_t ion_vel_z;
    real_t ion_temp;  // ion temperature in eV
    real_t ion_density;  // ion density in m^-3
    real_t sheath_potential = 3.1;
    real_t franck_condon_temp =
        2.0;    // temperature of dissociated molecules in eV
    real_t Rt;  // fraction of ions that are reflected
    real_t Rr;  // fraction of reflected particles that are backscattered
};

}  // namespace neptune

#endif