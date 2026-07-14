#ifndef NEPTUNE_SIMPLE_COLLISIONS_HPP
#define NEPTUNE_SIMPLE_COLLISIONS_HPP

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/physics/physicalconstants.hpp"

//
// This file contains collision processes with simple models
//
// rates from
// N. Horsten, W. Dekeyser, G. Samaey, and M. Baelmans,
// Comparison of fluid neutral models for one-dimensional plasma edge
// modeling with a finite volume solution of the Boltzmann equation
// (who got them from)
// R. J. Goldston and P. H. Rutherford, Introduction to Plasma
// Physics (CRC Press, 1995)
//

namespace neptune {

class SimpleChargeExchangeCollision : public ChargeExchangeCollision {
   public:
    SimpleChargeExchangeCollision(Species background_species,
                                  Species particle_species)
        : background_species_(background_species),
          particle_species_(particle_species)
    {}

    Species background_species() const override { return background_species_; }

    Species particle_species() const override { return particle_species_; }

    real_t collision_rate(const Particle& p, const Background& b) const override
    {
        real_t b_dens = b.density(background_species(), p.geometry.cell_id);
        real_t b_temp = b.temperature(background_species(), p.geometry.cell_id);
        return b_dens * 3.2e-15 * std::sqrt(b_temp / 0.026);
    }

    Vec3 sample_post_collision_velocity(Particle& p,
                                        const Background& b) const override
    {
        real_t vel_stddev = std::sqrt(this->post_collision_velocity_var(p, b));

        Vec3 b_vel = b.velocity(background_species(), p.geometry.cell_id);

        std::normal_distribution<real_t> dist_x(b_vel.x, vel_stddev);
        std::normal_distribution<real_t> dist_y(b_vel.y, vel_stddev);
        std::normal_distribution<real_t> dist_z(b_vel.z, vel_stddev);

        return Vec3(dist_x(p.rng), dist_y(p.rng), dist_z(p.rng));
    }

    Vec3 post_collision_velocity_mean(const Particle& p,
                                      const Background& b) const override
    {
        return b.velocity(background_species(), p.geometry.cell_id);
    }

    real_t post_collision_velocity_var(const Particle& p,
                                       const Background& b) const override
    {
        real_t b_temp = b.temperature(background_species(), p.geometry.cell_id);
        return b_temp * constants::eV_to_Joules /
               (constants::m_p *
                SPECIES_MASS[static_cast<index_t>(background_species_)]);
    }

   private:
    Species background_species_;
    Species particle_species_;
};

class SimpleIonisationCollision : public IonisationCollision {
   public:
    SimpleIonisationCollision(Species background_species,
                              Species particle_species)
        : background_species_(background_species),
          particle_species_(particle_species)
    {}

    Species background_species() const override { return background_species_; }

    Species particle_species() const override { return particle_species_; }

    real_t collision_rate(const Particle& p, const Background& b) const override
    {
        // K_i = 2.0 * 1e-13 / (6.0 + te / (13.6))\
        //     * np.sqrt(te / (13.6))\
        //     * np.exp(-13.6 / (te))
        // R_i = ne * K_i

        real_t e_dens = b.density(Species::e, p.geometry.cell_id);
        real_t e_temp = b.temperature(Species::e, p.geometry.cell_id);
        return e_dens * 2.0 * 1e-13 / (6.0 + e_temp / (13.6)) *
               std::sqrt(e_temp / (13.6)) * std::exp(-13.6 / (e_temp));
    }

   private:
    Species background_species_;
    Species particle_species_;
};

}  // namespace neptune

#endif