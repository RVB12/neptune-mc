#ifndef NEPTUNE_AMJUEL_COLLISIONS_HPP
#define NEPTUNE_AMJUEL_COLLISIONS_HPP

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/physics/physicalconstants.hpp"

//
// This file contains collision processes with rates from the AMJUEL database
//

namespace neptune {

/**
 * This is a charge exchange collision with rates from the AMJUEL database.
 * You should provide a 9x9 matrix of coefficients that define the collision rate
 * as a function of the logarithm of the particle energy and the logarithm of the
 * background temperature. i.e. H.3 in the AMJUEL database.
 */
class AmjuelChargeExchangeCollision : public ChargeExchangeCollision {
   public:
    AmjuelChargeExchangeCollision(Species background_species,
                                  Species particle_species,
                                  const vector<real_t>& coef_matrix)
        : background_species_(background_species),
          particle_species_(particle_species)
    {
        for (int i = 0; i < 81; ++i) {
            coef_matrix_[i] = coef_matrix[i];
        }
    }

    Species background_species() const override { return background_species_; }

    Species particle_species() const override { return particle_species_; }

    real_t collision_rate(const Particle& p, const Background& b) const override
    {
        Vec3 b_vel = b.velocity(background_species(), p.geometry.cell_id);
        Vec3 rel_vel = p.geometry.velocity - b_vel;

        real_t p_energy =
            0.5 * SPECIES_MASS[static_cast<index_t>(particle_species_)] *
            constants::m_p * rel_vel.norm2() / constants::eV_to_Joules;
        real_t b_temp = b.temperature(background_species(), p.geometry.cell_id);
        real_t b_dens = b.density(background_species(), p.geometry.cell_id);

        if(p_energy <= 0.0 or b_temp <= 0.0 or b_dens <= 0.0) {
            return 0.0;
        }

        real_t logT = std::log(p_energy);
        real_t logE = std::log(b_temp);

        // evaluate sum_i sum_j coef_matrix_[i][j] * logT^i * logE^j
        real_t rate = 0.0;
        for (int i = 8; i >= 0; --i) {
            real_t row = 0.0;
            for (int j = 8; j >= 0; --j) {
                row = row * logE + coef_matrix_[i * 9 + j];
            }
            rate = rate * logT + row;
        }

        // rate is in cm^3/s so also multiply by 1e-6 to convert to m^3/s
        real_t sigma_v = std::exp(rate) * 1e-6;

        return b_dens * sigma_v;
    }

    // Sampling the post-collisional velocity is still just a sample from a
    // Maxwellian, same as for the simple charge exchange collision
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
    real_t coef_matrix_[81];
};


/**
 * This is a charge exchange collision with rates from the AMJUEL database.
 * You should provide a 9x9 matrix of coefficients that define the collision rate
 * as a function of the logarithm of the electron density and the logarithm of the
 * electron temperature. i.e. H.4 in the AMJUEL database.
 */
class AmjuelIonisationCollision : public IonisationCollision {
   public:
    AmjuelIonisationCollision(Species background_species,
                              Species particle_species,
                              const vector<real_t>& coef_matrix)
        : background_species_(background_species),
          particle_species_(particle_species)
    {
        for (int i = 0; i < 81; ++i) {
            coef_matrix_[i] = coef_matrix[i];
        }
    }

    Species background_species() const override { return background_species_; }

    Species particle_species() const override { return particle_species_; }

    real_t collision_rate(const Particle& p, const Background& b) const override
    {
        real_t b_temp = b.temperature(Species::e, p.geometry.cell_id);
        real_t b_dens = b.density(Species::e, p.geometry.cell_id);

        if(b_temp <= 0.0 or b_dens <= 0.0) {
            return 0.0;
        }

        real_t logT = std::log(b_temp);
        real_t logE = std::log(b_dens/1e8/1e6); // convert to cm^-3 and then to 10^8 cm^-3

        // evaluate sum_i sum_j coef_matrix_[i][j] * logT^i * logE^j
        real_t rate = 0.0;
        for (int i = 8; i >= 0; --i) {
            real_t row = 0.0;
            for (int j = 8; j >= 0; --j) {
                row = row * logE + coef_matrix_[i * 9 + j];
            }
            rate = rate * logT + row;
        }

        // rate is in cm^3/s so also multiply by 1e-6 to convert to m^3/s
        real_t sigma_v = std::exp(rate) * 1e-6;

        return b_dens * sigma_v;
    }

   private:
    Species background_species_;
    Species particle_species_;
    real_t coef_matrix_[81];
};

}  // namespace neptune

#endif