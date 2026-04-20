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

class SimpleChargeExchangeCollision_D : public ChargeExchangeCollision {
   public:
    SimpleChargeExchangeCollision_D() = default;

    Species background_species() const override
    {
        return Species::Dplus;
    }

    Species particle_species() const override
    {
        return Species::D;
    }

    real_t collision_rate(real_t b_dens,
                          Vec3 b_vel,
                          real_t b_temp,
                          const Particle& p) const override
    {
        return b_dens * 3.2e-15 * std::sqrt(b_temp / 0.026);
    }

    Vec3 sample_post_collision_velocity(real_t b_dens,
                                        Vec3 b_vel,
                                        real_t b_temp,
                                        Particle& p) const override
    {
        real_t vel_stddev =
            std::sqrt(this->post_collision_velocity_var(b_dens, b_vel, b_temp, p));
        std::normal_distribution<real_t> dist_x(b_vel.x, vel_stddev);
        std::normal_distribution<real_t> dist_y(b_vel.y, vel_stddev);
        std::normal_distribution<real_t> dist_z(b_vel.z, vel_stddev);

        return Vec3(dist_x(p.rng), dist_y(p.rng), dist_z(p.rng));
    }

    Vec3 post_collision_velocity_mean(real_t b_dens,
                                      Vec3 b_vel,
                                      real_t b_temp,
                                      const Particle& p) const override
    {
        return b_vel;
    }

    real_t post_collision_velocity_var(real_t b_dens,
                                       Vec3 b_vel,
                                       real_t b_temp,
                                       const Particle& p) const override
    {
        return b_temp * constants::eV_to_Joules /
               (constants::m_p *
                SPECIES_MASS[static_cast<index_t>(Species::Dplus)]);
    }
};


class SimpleIonisationCollision_D : public IonisationCollision {
   public:
    SimpleIonisationCollision_D() = default;

    Species background_species() const override
    {
        return Species::Dplus;
    }

    Species particle_species() const override
    {
        return Species::D;
    }

    real_t collision_rate(real_t b_dens,
                          Vec3 b_vel,
                          real_t b_temp,
                          const Particle& p) const override
    {
        return b_dens * 2.0e-13/(6 + b_temp/13.6) * std::sqrt(b_temp / 13.6) * std::exp(-13.6/b_temp);
    }
};

}  // namespace neptune

#endif