#ifndef EIRON2_LINESOURCE_HPP
#define EIRON2_LINESOURCE_HPP

#include <cassert>
#include <random>

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vec2.hpp"
#include "eiron2/physics/particle.hpp"
#include "eiron2/physics/source.hpp"

namespace eiron2 {

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
    LineSource()
        : Source()
    {}

    LineSource(index_t max_particles_, Vec2 pos_a_, Vec2 pos_b_)
        : Source(max_particles_), pos_a(pos_a_), pos_b(pos_b_)
    {}

    LineSource(index_t max_particles_,
               Vec2 pos_a_,
               Vec2 pos_b_,
               real_t angle_min_,
               real_t angle_max_)
        : Source(max_particles_),
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
    Particle generate_particle();

    Vec2 pos_a;                // start of the line that defines this source
    Vec2 pos_b;                // end of the line that defines this source
    real_t angle_min = -M_PI/2;  // minimum angle of the velocity vector
    real_t angle_max = M_PI/2;   // maximum angle of the velocity vector
};

}  // namespace eiron2

#endif