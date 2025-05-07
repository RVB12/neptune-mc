#ifndef NEPTUNE_LINESOURCE_HPP
#define NEPTUNE_LINESOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/source.hpp"

namespace neptune {

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
class LineSource : public Source {
   public:
    LineSource(index_t num_particles_,
               real_t intensity_,
               Vec3 pos_a_,
               Vec3 pos_b_,
               Vec3 vel_mean_,
               real_t vel_stddev_)
        : Source(num_particles_, intensity_), pos_a(pos_a_), pos_b(pos_b_), vel_mean(vel_mean_), vel_stddev(vel_stddev_)
    {}

    /**
     * @brief Generate a particle
     *
     * Note: This will not set the cell_idx!! You need to do this afterwards!
     */
    Particle generate_particle(index_t i_particle);

    Vec3 pos_a;                    // start of the line that defines this source
    Vec3 pos_b;                    // end of the line that defines this source
    Vec3 vel_mean; // Mean velocity of the particles
    real_t vel_stddev; // Standard deviation of the velocity of the particles
};

}  // namespace neptune

#endif