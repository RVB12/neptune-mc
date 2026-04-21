#include "neptune/physics/pointsource.hpp"
#include "neptune/physics/physicalconstants.hpp"

#include <iostream>

namespace neptune {

Particle PointSource::generate_particle(index_t i_particle, const Mesh& mesh, const Background& background)
{
    if (i_particle >= num_particles()) {
        throw std::runtime_error(
            "PointSource::generate_particle: "
            "Tried to generate more particles than "
            "allowed by num_particles");
    }
    Particle p(i_particle + particle_id_offset_);

    if(cell_idx == -1){
        cell_idx = mesh.get_cell_index(pos); // Find cell_idx if it was not yet set
        if(cell_idx == -1) {
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