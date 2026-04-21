#include "neptune/estimation/estimators.hpp"

#include "neptune/finitevolume/advectiondiffusionsolver.hpp"
#include "neptune/physics/physicalconstants.hpp"
namespace neptune {

void Estimator_old::add_track_length(const Particle& particle,
                                 const Particle& particle_old,
                                 const Mesh& mesh,
                                 const real_t R_i)
{
    real_t weight_old = particle_old.weight;
    real_t t_old = particle_old.t;
    index_t cell_id_old = particle_old.geometry.cell_id;
    real_t volume = mesh.get_volume(cell_id_old);

    real_t dt = particle.t - t_old;
    real_t s;
    if(R_i * dt < 1.0e-4) {
        // Avoid division by zero if ionisation rate is very low
        s = weight_old * dt / volume * particle.intensity;
    } else {
        s = weight_old * (1 - std::exp(-R_i * dt))/R_i / volume * particle.intensity;
    }

    density[cell_id_old] += s;

    momentum_x[cell_id_old] += particle_old.geometry.velocity.x * s;
    momentum_y[cell_id_old] += particle_old.geometry.velocity.y * s;
    momentum_z[cell_id_old] += particle_old.geometry.velocity.z * s;

    energy[cell_id_old] +=
        s * (0.5 * SPECIES_MASS[static_cast<size_t>(particle_old.species)] * particle_old.geometry.velocity.norm2()) *
        constants::m_p / constants::eV_to_Joules;

    energy_x[cell_id_old] +=
        s *
        (0.5 * SPECIES_MASS[static_cast<size_t>(particle_old.species)] * particle_old.geometry.velocity.x *
         particle_old.geometry.velocity.x) *
        constants::m_p / constants::eV_to_Joules;

    moment2_x[cell_id_old] +=
    s * particle_old.geometry.velocity.x * particle_old.geometry.velocity.x;

    moment3_x[cell_id_old] +=
    s * particle_old.geometry.velocity.x * particle_old.geometry.velocity.x * particle_old.geometry.velocity.x;
}

void Estimator_old::add_final(const Particle& particle, const Mesh& mesh)
{
    index_t cell_id = particle.geometry.cell_id;
    real_t volume = mesh.get_volume(cell_id);
    density_final[cell_id] += particle.weight / volume * particle.intensity;

    momentum_x_final[cell_id] += particle.geometry.velocity.x *
                                 particle.weight / volume * particle.intensity;
    momentum_y_final[cell_id] += particle.geometry.velocity.y *
                                 particle.weight / volume * particle.intensity;
    momentum_z_final[cell_id] += particle.geometry.velocity.z *
                                 particle.weight / volume * particle.intensity;
    energy_final[cell_id] +=
        (0.5 * SPECIES_MASS[static_cast<size_t>(particle.species)] * particle.geometry.velocity.norm2()) *
        particle.weight / volume * particle.intensity;

    energy_x_final[cell_id] +=
        (0.5 * SPECIES_MASS[static_cast<size_t>(particle.species)] * particle.geometry.velocity.x *
         particle.geometry.velocity.x) *
        particle.weight / volume * particle.intensity;
}

void Estimator_old::combine(const Estimator_old& other)
{
    for (index_t i = 0; i < density.size(); i++) {
        density[i] += other.density[i];
        momentum_x[i] += other.momentum_x[i];
        momentum_y[i] += other.momentum_y[i];
        momentum_z[i] += other.momentum_z[i];
        energy[i] += other.energy[i];
        energy_x[i] += other.energy_x[i];
        moment2_x[i] += other.moment2_x[i];
        moment3_x[i] += other.moment3_x[i];
        density_final[i] += other.density_final[i];
        momentum_x_final[i] += other.momentum_x_final[i];
        momentum_y_final[i] += other.momentum_y_final[i];
        momentum_z_final[i] += other.momentum_z_final[i];
        energy_final[i] += other.energy_final[i];
        energy_x_final[i] += other.energy_x_final[i];
    }
}

void Estimator_old::combine(const FluidEstimator& fluid_estimator)
{
    for (index_t i = 0; i < density.size(); i++) {
        density[i] += fluid_estimator.fluid_density_integral[i];
        momentum_x[i] += fluid_estimator.fluid_momentum_x_integral[i];
        momentum_y[i] += fluid_estimator.fluid_momentum_y_integral[i];
        momentum_z[i] += fluid_estimator.fluid_momentum_z_integral[i];
        energy[i] += fluid_estimator.fluid_energy_integral[i];
        energy_x[i] += fluid_estimator.fluid_energy_x_integral[i];
    }
}

void Estimator_old::combine(const ArcmcEstimator& other)
{
    for (index_t i = 0; i < density.size(); i++) {
        density[i] += other.density[i];
        momentum_x[i] += other.momentum_x[i];
        momentum_y[i] += other.momentum_y[i];
        momentum_z[i] += other.momentum_z[i];
        energy[i] += other.energy[i];
        energy_x[i] += other.energy_x[i];
        density_final[i] += other.density_final[i];
        momentum_x_final[i] += other.momentum_x_final[i];
        momentum_y_final[i] += other.momentum_y_final[i];
        momentum_z_final[i] += other.momentum_z_final[i];
        energy_final[i] += other.energy_final[i];
        energy_x_final[i] += other.energy_x_final[i];
    }
}

void FluidEstimator::add_fluid_step(const Particle& particle,
                                    const Particle& particle_old,
                                    const Mesh& mesh)
{
    real_t weight_old = particle_old.weight;
    real_t t_old = particle_old.t;
    index_t cell_id_old = particle_old.geometry.cell_id;
    real_t volume = mesh.get_volume(cell_id_old);

    real_t s = weight_old / volume * particle.intensity;
    initial_density[cell_id_old] += s;

    index_t cell_id = particle.geometry.cell_id;
    real_t volume_new = mesh.get_volume(cell_id);
    real_t s2 = particle.weight / volume_new * particle.intensity;
    end_density[cell_id] += s2;

    dt += particle.t - t_old;
    n_particles++;
}

void FluidEstimator::combine(const FluidEstimator& other)
{
    for (index_t i = 0; i < initial_density.size(); i++) {
        initial_density[i] += other.initial_density[i];
        end_density[i] += other.end_density[i];
    }
    dt += other.dt;
    n_particles += other.n_particles;
}

void FluidEstimator::finalize(
    Mesh& mesh,
    const vector<std::unique_ptr<FluidBoundaryCondition>>& boundary_conditions,
    const vector<real_t>& advection_x,
    const vector<real_t>& advection_y,
    const vector<real_t>& advection_z,
    const vector<real_t>& diffusion_rate,
    const vector<real_t>& diffusion_factor,
    const vector<real_t>& ionisation_rate)
{
    // Determine the average timestep
    dt /= n_particles;

    fluid_density_integral.resize(initial_density.size());
    fluid_momentum_x_integral.resize(initial_density.size());
    fluid_momentum_y_integral.resize(initial_density.size());
    fluid_momentum_z_integral.resize(initial_density.size());
    fluid_energy_integral.resize(initial_density.size());
    fluid_energy_x_integral.resize(initial_density.size());

    real_t dt_fluid2 = dt;

    std::cout << "Running fluid solver with dt = " << dt_fluid2 << std::endl;

    advection_diffusion_solver(
        mesh, boundary_conditions, advection_x, advection_y, advection_z,
        diffusion_rate, diffusion_factor, ionisation_rate, dt, dt_fluid2, initial_density,
        fluid_density_integral, fluid_momentum_x_integral,
        fluid_momentum_y_integral, fluid_momentum_z_integral,
        fluid_energy_integral, fluid_energy_x_integral);
}

// This track length estimator is almost identical to the normal estimator
// except that we correct the energy estimate by a factor of 1/W to account for
// difference between the effective velocity and the actual velocity.
void ArcmcEstimator::add_track_length(const Particle& particle,
                                      const Particle& particle_old,
                                      const Mesh& mesh,
                                      real_t R_i,
                                      real_t V2_eff,
                                      real_t V2_eff_x)
{
    real_t weight_old = particle_old.weight;
    real_t t_old = particle_old.t;
    index_t cell_id_old = particle_old.geometry.cell_id;
    real_t volume = mesh.get_volume(cell_id_old);

    real_t dt = particle.t - t_old;
    real_t s;
    if(R_i * dt < 1.0e-4) {
        // Avoid division by zero if ionisation rate is very low
        s = weight_old * dt / volume * particle.intensity;
    } else {
        s = weight_old * (1 - std::exp(-R_i * dt))/R_i / volume * particle.intensity;
    }

    density[cell_id_old] += s;

    momentum_x[cell_id_old] += particle_old.geometry.velocity.x * s;
    momentum_y[cell_id_old] += particle_old.geometry.velocity.y * s;
    momentum_z[cell_id_old] += particle_old.geometry.velocity.z * s;

    energy[cell_id_old] += s * (0.5 * SPECIES_MASS[static_cast<size_t>(particle_old.species)] * V2_eff) *
                           constants::m_p / constants::eV_to_Joules;

    energy_x[cell_id_old] += s * (0.5 * SPECIES_MASS[static_cast<size_t>(particle_old.species)] * V2_eff_x) *
                             constants::m_p / constants::eV_to_Joules;
}

void ArcmcEstimator::add_final(const Particle& particle,
                               const Mesh& mesh,
                               real_t V2_eff,
                               real_t V2_eff_x)
{
    index_t cell_id = particle.geometry.cell_id;
    real_t volume = mesh.get_volume(cell_id);
    density_final[cell_id] += particle.weight / volume * particle.intensity;

    momentum_x_final[cell_id] += particle.geometry.velocity.x *
                                 particle.weight / volume * particle.intensity;
    momentum_y_final[cell_id] += particle.geometry.velocity.y *
                                 particle.weight / volume * particle.intensity;
    momentum_z_final[cell_id] += particle.geometry.velocity.z *
                                 particle.weight / volume * particle.intensity;
    energy_final[cell_id] += (0.5 * SPECIES_MASS[static_cast<size_t>(particle.species)] * V2_eff) * particle.weight /
                             volume * particle.intensity;

    energy_x_final[cell_id] += (0.5 * SPECIES_MASS[static_cast<size_t>(particle.species)] * V2_eff_x) *
                               particle.weight / volume * particle.intensity;
}

void ArcmcEstimator::combine(const ArcmcEstimator& other)
{
    for (index_t i = 0; i < density.size(); i++) {
        density[i] += other.density[i];
        momentum_x[i] += other.momentum_x[i];
        momentum_y[i] += other.momentum_y[i];
        momentum_z[i] += other.momentum_z[i];
        energy[i] += other.energy[i];
        energy_x[i] += other.energy_x[i];
        density_final[i] += other.density_final[i];
        momentum_x_final[i] += other.momentum_x_final[i];
        momentum_y_final[i] += other.momentum_y_final[i];
        momentum_z_final[i] += other.momentum_z_final[i];
        energy_final[i] += other.energy_final[i];
        energy_x_final[i] += other.energy_x_final[i];
    }
}

}  // namespace neptune