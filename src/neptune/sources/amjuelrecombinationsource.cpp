#include "neptune/sources/amjuelrecombinationsource.hpp"

#include <iostream>

#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

AmjuelRecombinationSource::AmjuelRecombinationSource(index_t num_particles,
                                         index_t cell_idx_,
                                         Species background_species_,
                                         Species species_,
                                         const Mesh& mesh,
                                         const Background& background,
                                         const vector<real_t>& coef_matrix)
    : Source(num_particles),
      cell_idx(cell_idx_),
      background_species(background_species_),
      species(species_)
{
    real_t temperature = background.temperature(Species::e, cell_idx);
    real_t n_e = background.density(Species::e, cell_idx);


    real_t logT = std::log(temperature);
    real_t logE = std::log(n_e/1e8/1e6); // convert to cm^-3 and then to 10^8 cm^-3

    // evaluate sum_i sum_j coef_matrix_[i][j] * logT^i * logE^j
    real_t rate = 0.0;
    for (int i = 8; i >= 0; --i) {
        real_t row = 0.0;
        for (int j = 8; j >= 0; --j) {
            row = row * logE + coef_matrix[i * 9 + j];
        }
        rate = rate * logT + row;
    }

    // rate is in cm^3/s so also multiply by 1e-6 to convert to m^3/s
    real_t K_r = std::exp(rate) * 1e-6;
    real_t n_Dplus = background.density(background_species, cell_idx);
    real_t vol = mesh.get_volume(cell_idx);
    intensity_ = K_r * n_Dplus * n_e * vol;
}

Particle AmjuelRecombinationSource::generate_particle(uint64_t particle_id,
                                                uint64_t seed,
                                                const Mesh& mesh,
                                                const Background& background)
{
    Particle p(particle_id, seed, species);

    p.geometry.cell_id = cell_idx;
    p.geometry.boundary_id = -1;

    p.geometry.position = mesh.generate_random_position(cell_idx, p.rng);

    Vec3 vel_mean = background.velocity(background_species, cell_idx);
    real_t temperature = background.temperature(background_species, cell_idx);

    real_t vel_stddev = std::sqrt(
        constants::eV_to_Joules * temperature /
        (constants::m_p * SPECIES_MASS[static_cast<int>(background_species)]));

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