#include "neptune/physics/recombinationsource.hpp"

#include <iostream>

#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

RecombinationSource::RecombinationSource(index_t num_particles,
                                         index_t cell_idx_,
                                         const Mesh& mesh,
                                         const Background& background)
    : Source(num_particles), cell_idx(cell_idx_)
{
    real_t temperature = background.temperature(Species::Dplus, cell_idx);
    real_t K_r =
        0.7e-19 *
        std::sqrt(13.6 / temperature);  // Recombination rate coefficient
    real_t n_Dplus = background.density(Species::Dplus, cell_idx);
    real_t n_e = n_Dplus;  // Assuming quasi-neutrality, the electron density is
                           // equal to the D+ density
    real_t vol = mesh.get_volume(cell_idx);
    intensity_ = K_r * n_Dplus * n_e * vol;
}

Particle RecombinationSource::generate_particle(uint64_t particle_id,
                                                 rng_key_t key,
                                                 const Mesh& mesh,
                                                 const Background& background)
{
    Particle p(particle_id, key);

    p.geometry.cell_id = cell_idx;
    p.geometry.boundary_id = -1;

    p.geometry.position = mesh.generate_random_position(cell_idx, p.rng);

    Vec3 vel_mean = background.velocity(Species::Dplus, cell_idx);
    real_t temperature = background.temperature(Species::Dplus, cell_idx);

    real_t vel_stddev =
        std::sqrt(constants::eV_to_Joules * temperature / (constants::m_p *
                  SPECIES_MASS[static_cast<int>(Species::Dplus)]));

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