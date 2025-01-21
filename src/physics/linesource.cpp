#include "eiron2/physics/linesource.hpp"

#include <iostream>

namespace eiron2 {

Particle LineSource::generate_particle(index_t i_particle)
{
    if (i_particle >= num_particles()) {
        throw std::runtime_error(
            "LineSource::generate_particle: "
            "Tried to generate more particles than "
            "allowed by num_particles");
    }
    Particle p(i_particle + particle_id_offset_);

    // The position is determined as pos_a + u * (pos_b - pos_a),
    // where u is a uniform random number between 0 and 1.
    // This amounts to a random point on the line segment defined by pos_a
    // and pos_b.
    std::uniform_real_distribution<real_t> unif_dist(0., 1.);
    real_t u = unif_dist(p.rng);
    p.geometry.position = pos_a + u * (pos_b - pos_a);

    // The velocity is determined in two steps:
    // first, we generate a random angle between angle_min and angle_max
    // next, we draw the norm of the velocity from the maxwellian
    // distribution
    std::uniform_real_distribution<real_t> angle_dist(angle_min, angle_max);
    // The angle of the normal vector of the line segment
    real_t normal_angle = std::atan2(pos_b.y - pos_a.y, pos_b.x - pos_a.x) - M_PI/2;
    real_t angle = angle_dist(p.rng) + normal_angle;
    std::chi_squared_distribution<real_t> chi2_dist(3.);
    // todo: scale by sqrt(T/m)
    real_t v_nrm = std::sqrt(chi2_dist(p.rng));
 
    p.geometry.velocity.x = v_nrm * std::cos(angle);
    p.geometry.velocity.y = v_nrm * std::sin(angle);

    return p;
}

}  // namespace eiron2