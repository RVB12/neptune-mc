#include "neptune/sources/volumesource.hpp"

#include <iostream>

#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

Particle VolumeSource::generate_particle(
    uint64_t particle_id,
    uint64_t seed,
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions)
{
    Particle p(particle_id, seed, species);

    p.geometry.cell_id = cell_idx;
    p.geometry.boundary_id = -1;

    p.geometry.position = mesh.generate_random_position(cell_idx, p.rng);

    std::normal_distribution<real_t> dist_vx(vel_mean.x, vel_stddev);
    std::normal_distribution<real_t> dist_vy(vel_mean.y, vel_stddev);
    std::normal_distribution<real_t> dist_vz(vel_mean.z, vel_stddev);

    p.geometry.velocity = Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

    std::exponential_distribution<real_t> dist_t(1.0);
    p.u = dist_t(p.rng);

    p.intensity = intensity() / num_particles();

    return p;
}

}  // namespace neptune