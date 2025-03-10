#ifndef EIGRENE_POINTSOURCE_HPP
#define EIGRENE_POINTSOURCE_HPP

#include <cassert>
#include <random>

#include "eigrene/common/utils.hpp"
#include "eigrene/geometry/vec2.hpp"
#include "eigrene/physics/particle.hpp"
#include "eigrene/physics/source.hpp"

namespace eigrene {

/**
 * @brief A line source of particles
 *
 * The particles will be generated somewhere on a line segment defined by two
 * points. The velocity is drawn from a maxwellian distribution. It is possible
 * to restrict the possible angles of the velocity vector.
 *
 * By default, the velocity will only be to the right of the line segment:
 *
 */
class PointSource : public Source {
   public:
    PointSource(index_t num_particles_,
               real_t intensity_,
               Vec3 pos_,
               Vec3 vel_mean_,
               real_t vel_stddev_)
        : Source(num_particles_, intensity_), pos(pos_), vel_mean(vel_mean_), vel_stddev(vel_stddev_)
    {}

    /**
     * @brief Generate a particle
     *
     * Note: This will not set the cell_idx!! You need to do this afterwards!
     */
    Particle generate_particle(index_t i_particle);

    Vec3 pos; // Initial position of the particles
    Vec3 vel_mean; // Mean velocity of the particles
    real_t vel_stddev; // Standard deviation of the velocity of the particles
};

}  // namespace eigrene

#endif