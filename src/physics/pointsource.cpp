#include "neptune/physics/pointsource.hpp"

#include <iostream>

#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

Particle PointSource::generate_particle(uint64_t particle_id,
                                        rng_key_t key,
                                        const Mesh& mesh,
                                        const Background& background)
{
    Particle p(particle_id, key);

    if (cell_idx == -1) {
        cell_idx =
            mesh.get_cell_index(pos);  // Find cell_idx if it was not yet set
        if (cell_idx == -1) {
            throw std::runtime_error(
                "PointSource::generate_particle: "
                "Could not find cell index for position");
        }
    }

    p.geometry.cell_id = cell_idx;
    p.geometry.position = pos;

    std::normal_distribution<real_t> dist_vx(vel_mean.x, vel_stddev);
    std::normal_distribution<real_t> dist_vy(vel_mean.y, vel_stddev);
    std::normal_distribution<real_t> dist_vz(vel_mean.z, vel_stddev);

    p.geometry.velocity = Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

    std::exponential_distribution<real_t> dist_t(1.0);
    p.u = dist_t(p.rng);

    return p;
}

}  // namespace neptune