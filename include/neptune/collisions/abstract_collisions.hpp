#ifndef NEPTUNE_ABSTRACTCOLLISIONS_HPP
#define NEPTUNE_ABSTRACTCOLLISIONS_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/vec3.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/species.hpp"

namespace neptune {

class Collision {
    public:
    // Species of the background particle involved in the collision
    virtual Species background_species() const = 0;

    // Species of the particle involved in the collision
    virtual Species particle_species() const = 0;

    /**
     * Calculate the collision rate (in s^-1)
     *
     * @param b_dens The background density of the species involved in the
     * collision (m^-3)
     * @param b_vel The background velocity of the species involved in the
     * collision (m/s)
     * @param b_temp The background temperature of the species involved
     * in the collision (eV)
     * @param p The particle for which to calculate the collision rate
     * @return The collision rate in s^-1
     */
    virtual real_t collision_rate(real_t b_dens,
                                  Vec3 b_vel,
                                  real_t b_temp,
                                  const Particle& p) const = 0;
};

class ChargeExchangeCollision: public Collision {
   public:

    /**
     * Calculate the mean post-collision velocity of a particle
     *
     * @param b_dens The background density of the species involved in the
     * collision (m^-3)
     * @param b_vel The background velocity of the species involved in the
     * collision (m/s)
     * @param b_temp The background temperature of the species involved
     * in the collision (eV)
     * @param p The particle for which to calculate the mean post-collision
     * velocity
     * @return The mean post-collision velocity of the particle (m/s)
     */
    virtual Vec3 post_collision_velocity_mean(real_t b_dens,
                                              Vec3 b_vel,
                                              real_t b_temp,
                                              const Particle& p) const = 0;
    /**
     * Calculate the variance of the post-collision velocity of a particle
     *
     * @param b_dens The background density of the species involved in the
     * collision (m^-3)
     * @param b_vel The background velocity of the species involved in the
     * collision (m/s)
     * @param b_temp The background temperature of the species involved
     * in the collision (eV)
     * @param p The particle for which to calculate the variance post-collision
     * velocity
     * @return The variance post-collision velocity of the particle (m/s)
     */
    virtual real_t post_collision_velocity_var(real_t b_dens,
                                               Vec3 b_vel,
                                               real_t b_temp,
                                               const Particle& p) const = 0;

    /**
     * Sample the post-collision velocity of a particle
     *
     * @param b_dens The background density of the species involved in the
     * collision (m^-3)
     * @param b_vel The background velocity of the species involved in the
     * collision (m/s)
     * @param b_temp The background temperature of the species involved
     * in the collision (eV)
     * @param p The particle for which to sample the post-collision velocity
     * @return The post-collision velocity of the particle (m/s)
     */
    virtual Vec3 sample_post_collision_velocity(real_t b_dens,
                                                Vec3 b_vel,
                                                real_t b_temp,
                                                Particle& p) const = 0;
};

class IonisationCollision: public Collision {
    // Ionisation has nothing to add to the base Collision class,
    // but we define it as a separate class for clarity
};

}  // namespace neptune

#endif