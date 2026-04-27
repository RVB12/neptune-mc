#ifndef NEPTUNE_ABSTRACTCOLLISIONS_HPP
#define NEPTUNE_ABSTRACTCOLLISIONS_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/vec3.hpp"
#include "neptune/physics/background.hpp"
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
     * @param p The particle for which to calculate the collision rate
     * @param b The background for which to calculate the collision rate
     * @return The collision rate in s^-1
     */
    virtual real_t collision_rate(const Particle& p,
                                  const Background& b) const = 0;
};

class ChargeExchangeCollision : public Collision {
   public:
    /**
     * Calculate the mean post-collision velocity of a particle
     *
     * @param p The particle for which to calculate the mean post-collision
     * velocity
     * @param b The background for which to calculate the post-collision
     * velocity
     * @return The mean post-collision velocity of the particle (m/s)
     */
    virtual Vec3 post_collision_velocity_mean(const Particle& p,
                                              const Background& b) const = 0;
    /**
     * Calculate the variance of the post-collision velocity of a particle
     *
     * @param p The particle for which to calculate the variance post-collision
     * velocity
     * @param b The background for which to calculate the variance
     * post-collision velocity
     * @return The variance post-collision velocity of the particle (m/s)
     */
    virtual real_t post_collision_velocity_var(const Particle& p,
                                               const Background& b) const = 0;

    /**
     * Sample the post-collision velocity of a particle
     *
     * @param p The particle for which to sample the post-collision velocity
     * @param b The background for which to sample the post-collision velocity
     * @return The post-collision velocity of the particle (m/s)
     */
    virtual Vec3 sample_post_collision_velocity(Particle& p,
                                                const Background& b) const = 0;
};

class IonisationCollision : public Collision {
    // Ionisation has nothing to add to the base Collision class,
    // but we define it as a separate class for clarity
};

}  // namespace neptune

#endif