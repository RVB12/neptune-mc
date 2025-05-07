#include "neptune/physics/linesource.hpp"
#include "neptune/physics/physicalconstants.hpp"

#include <iostream>

namespace neptune {

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

    std::normal_distribution<real_t> dist_vx(vel_mean.x, vel_stddev);
    std::normal_distribution<real_t> dist_vy(vel_mean.y, vel_stddev);
    std::normal_distribution<real_t> dist_vz(vel_mean.z, vel_stddev);

    p.geometry.velocity = Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

    return p;
}

}  // namespace neptune