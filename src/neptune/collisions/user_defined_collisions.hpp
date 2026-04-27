#ifndef NEPTUNE_USER_DEFINED_COLLISIONS_HPP
#define NEPTUNE_USER_DEFINED_COLLISIONS_HPP

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/physics/constants.hpp"

//
// This file contains collision processes with user defined rates and post-collision velocity distributions
//

namespace neptune {

class UserDefinedChargeExchangeCollision : public ChargeExchangeCollision {
   public:
    UserDefinedChargeExchangeCollision(
        Species background_species,
        Species particle_species,
        vector<real_t> collision_rate,
        vector<real_t> bgk_mean_x,
        vector<real_t> bgk_mean_y,
        vector<real_t> bgk_mean_z,
        vector<real_t> bgk_var)
        : background_species_(background_species),
          particle_species_(particle_species),
          collision_rate(collision_rate),
          bgk_mean_x(bgk_mean_x),
          bgk_mean_y(bgk_mean_y),
          bgk_mean_z(bgk_mean_z),
          bgk_var(bgk_var)
    {}

    Species background_species() const override
    {
        return background_species_;
    }

    Species particle_species() const override
    {
        return particle_species_;
    }

     real_t collision_rate(real_t b_dens,
                          Vec3 b_vel,
                          real_t b_temp,
                          const Particle& p) const override
    {
        return collision_rate[p.geometry.cell_id];
    }

    Vec3 sample_post_collision_velocity(real_t b_dens,
                                        Vec3 b_vel,
                                        real_t b_temp,
                                        Particle& p) const override
    {
        real_t vel_stddev =
            std::sqrt(bgk_var[p.geometry.cell_id]);
        Vec3 mean_vel(bgk_mean_x[p.geometry.cell_id],
                      bgk_mean_y[p.geometry.cell_id],
                      bgk_mean_z[p.geometry.cell_id]);
        std::normal_distribution<real_t> dist_x(mean_vel.x, vel_stddev);
        std::normal_distribution<real_t> dist_y(mean_vel.y, vel_stddev);
        std::normal_distribution<real_t> dist_z(mean_vel.z, vel_stddev);

        return Vec3(dist_x(p.rng), dist_y(p.rng), dist_z(p.rng));
    }

    Vec3 post_collision_velocity_mean(real_t b_dens,
                                      Vec3 b_vel,
                                      real_t b_temp,
                                      const Particle& p) const override
    {
        return Vec3(bgk_mean_x[p.geometry.cell_id],
                    bgk_mean_y[p.geometry.cell_id],
                    bgk_mean_z[p.geometry.cell_id]);
    }

    real_t post_collision_velocity_var(real_t b_dens,
                                       Vec3 b_vel,
                                       real_t b_temp,
                                       const Particle& p) const override
    {
        return bgk_var[p.geometry.cell_id];
    }

    private:
        Species background_species_;
        Species particle_species_;
        vector<real_t> collision_rate;
        vector<real_t> bgk_mean_x;
        vector<real_t> bgk_mean_y;
        vector<real_t> bgk_mean_z;
        vector<real_t> bgk_var;
};

class UserDefinedIonisationCollision : public IonisationCollision {
   public:
    UserDefinedIonisation(
        Species background_species,
        Species particle_species,
        vector<real_t> collision_rate)
        : background_species_(background_species),
          particle_species_(particle_species),
          collision_rate(collision_rate)
    {}

    Species background_species() const override
    {
        return background_species_;
    }

    Species particle_species() const override
    {
        return particle_species_;
    }

     real_t collision_rate(real_t b_dens,
                          Vec3 b_vel,
                          real_t b_temp,
                          const Particle& p) const override
    {
        return collision_rate[p.geometry.cell_id];
    }

    private:
        Species background_species_;
        Species particle_species_;
        vector<real_t> collision_rate;
};

}  // namespace neptune

#endif