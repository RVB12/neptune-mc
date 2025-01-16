#ifndef EIRON2_LINESOURCE_HPP
#define EIRON2_LINESOURCE_HPP

#include <cassert>
#include <random>

#include "eiron2/common/vec2.hpp"
#include "eiron2/physics/particle.hpp"
#include "eiron2/common/utils.hpp"

namespace eiron2 {

/**
 * @brief A line source of particles
 *
 * The particles will be generated somewhere on a line segment defined by two
 * points. The velocity is drawn from a maxwellian distribution. It is possible
 * to restrict the possible angles of the velocity vector.
 *
 */
class LineSource {
   public:
    LineSource(Vec2 pos_a_, Vec2 pos_b_, index_t num_particles_)
        : pos_a(pos_a_), pos_b(pos_b_), num_particles(num_particles_)
    {
        assert(num_particles > 0);
    }

    LineSource(Vec2 pos_a_,
               Vec2 pos_b_,
               index_t num_particles_,
               real_t angle_min_,
               real_t angle_max_)
        : pos_a(pos_a_),
          pos_b(pos_b_),
          angle_min(angle_min_),
          angle_max(angle_max_),
          num_particles(num_particles_)
    {
        assert(num_particles > 0);
    }

    /**
     * @brief Generate a particle
     *
     * Note: This will not set the cell_idx!! You need to do this afterwards!
     */
    template <std::uniform_random_bit_generator RNG>
    Particle generate_particle(RNG rng) const
    {
        Particle p;

        // The position is determined as pos_a + u * (pos_b - pos_a),
        // where u is a uniform random number between 0 and 1.
        // This amounts to a random point on the line segment defined by pos_a
        // and pos_b.
        std::uniform_real_distribution<real_t> unif_dist(0., 1.);
        real_t u = unif_dist(rng);
        p.geometry.position = pos_a + u * (pos_b - pos_a);

        // The velocity is determined in two steps:
        // first, we generate a random angle between angle_min and angle_max
        // next, we draw the norm of the velocity from the maxwellian
        // distribution
        std::uniform_real_distribution<real_t> angle_dist(angle_min, angle_max);
        real_t angle = angle_dist(rng);
        // The angle is currently relative to the line segment, but we want it
        // relative to the x-axis
        angle += std::atan2(pos_b.y - pos_a.y, pos_b.x - pos_a.x);
        std::chi_squared_distribution<real_t> chi2_dist(3, p.mass);
        real_t v_nrm = std::sqrt(chi2_dist(rng));

        p.geometry.velocity.x = v_nrm * std::cos(angle);
        p.geometry.velocity.y = v_nrm * std::sin(angle);

        return p;
    }

   private:
    Vec2 pos_a;                // start of the line that defines this source
    Vec2 pos_b;                // end of the line that defines this source
    real_t angle_min = -M_PI;  // minimum angle of the velocity vector
    real_t angle_max = M_PI;   // maximum angle of the velocity vector
    index_t num_particles;  // the number of particles to generate from this source
};

}  // namespace eiron2

#endif