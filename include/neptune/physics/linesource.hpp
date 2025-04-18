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
               Vec3 pos_b_)
        : Source(num_particles_, intensity_), pos_a(pos_a_), pos_b(pos_b_)
    {}

    LineSource(index_t num_particles_,
               real_t intensity_,
               Vec3 pos_a_,
               Vec3 pos_b_,
               real_t angle_min_,
               real_t angle_max_)
        : Source(num_particles_, intensity_),
          pos_a(pos_a_),
          pos_b(pos_b_),
          angle_min(angle_min_),
          angle_max(angle_max_)
    {}

    /**
     * @brief Generate a particle
     *
     * Note: This will not set the cell_idx!! You need to do this afterwards!
     */
    Particle generate_particle(index_t i_particle);

    Vec3 pos_a;                    // start of the line that defines this source
    Vec3 pos_b;                    // end of the line that defines this source
    real_t angle_min = -M_PI / 2;  // minimum angle of the velocity vector
    real_t angle_max = M_PI / 2;   // maximum angle of the velocity vector
};

}  // namespace neptune

#endif