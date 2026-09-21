#include "neptune/sources/pointsource.hpp"

#include <iostream>

#include "neptune/physics/physicalconstants.hpp"
#include "neptune/geometry/cart_to_tor.hpp"

namespace neptune {

Particle PointSource::generate_particle(
    uint64_t particle_id,
    uint64_t seed,
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions)
{
    Particle p(particle_id, seed, species);

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

    if(mesh.is_toroidal()){
        p.geometry.velocity = transform_tor_dir_to_cart(p.geometry.position, p.geometry.velocity);
    }

    std::exponential_distribution<real_t> dist_t(1.0);
    p.u = dist_t(p.rng);

    p.intensity = intensity() / num_particles();

    return p;
}

}  // namespace neptune