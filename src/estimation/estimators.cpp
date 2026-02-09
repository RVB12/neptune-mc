#include "neptune/estimation/estimators.hpp"

#include "neptune/finitevolume/advectiondiffusionsolver.hpp"
#include "neptune/physics/physicalconstants.hpp"
namespace neptune {

void Estimator::add_track_length(const Particle& particle,
                                 const Particle& particle_old,
                                 const Mesh& mesh)
{
    real_t weight_old = particle_old.weight;
    real_t t_old = particle_old.t;
    index_t cell_id_old = particle_old.geometry.cell_id;
    real_t volume = mesh.get_volume(cell_id_old);

    real_t s = (particle.t - t_old) * weight_old / volume * particle.intensity;

    density[cell_id_old] += s;

    momentum_x[cell_id_old] += particle_old.geometry.velocity.x * s;
    momentum_y[cell_id_old] += particle_old.geometry.velocity.y * s;
    momentum_z[cell_id_old] += particle_old.geometry.velocity.z * s;

    energy[cell_id_old] +=
        s * (0.5 * particle_old.mass * particle_old.geometry.velocity.norm2()) *
        constants::m_p / constants::eV_to_Joules;

    energy_x[cell_id_old] +=
        s * (0.5 * particle_old.mass * particle_old.geometry.velocity.x * particle_old.geometry.velocity.x) *
        constants::m_p / constants::eV_to_Joules;
}

void Estimator::add_final(const Particle& particle, const Mesh& mesh)
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
        (0.5 * particle.mass * particle.geometry.velocity.norm2()) *
        particle.weight / volume * particle.intensity;

    energy_x_final[cell_id] +=
        (0.5 * particle.mass * particle.geometry.velocity.x * particle.geometry.velocity.x) *
        particle.weight / volume * particle.intensity;
}

void Estimator::combine(const Estimator& other)
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

void Estimator::combine(const FluidEstimator& fluid_estimator)
{
    for (index_t i = 0; i < density.size(); i++) {
        density[i] += fluid_estimator.fluid_density_integral[i];
        // TODO: add momentum and energy contributions from fluid estimator
    }
}

void Estimator::combine(const ArcmcEstimator& other)
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
    density[cell_id_old] += s;

    dt += particle.t - t_old;
    n_particles++;
}

void FluidEstimator::combine(const FluidEstimator& other)
{
    for (index_t i = 0; i < density.size(); i++) {
        density[i] += other.density[i];
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
    const vector<real_t>& diffusion_factor)
{
    // Determine the average timestep
    dt /= n_particles;

    fluid_density_integral.resize(density.size());

    real_t dt_fluid2 = dt;

    advection_diffusion_solver(mesh, boundary_conditions, advection_x,
                               advection_y, advection_z, diffusion_rate,
                               diffusion_factor, dt, dt_fluid2, density,
                               fluid_density_integral);
}


// This track length estimator is almost identical to the normal estimator
// except that we correct the energy estimate by a factor of 1/W to account for difference
// between the effective velocity and the actual velocity.
void ArcmcEstimator::add_track_length(const Particle& particle,
                                 const Particle& particle_old,
                                 const Mesh& mesh,
                                 real_t W)
{
    real_t weight_old = particle_old.weight;
    real_t t_old = particle_old.t;
    index_t cell_id_old = particle_old.geometry.cell_id;
    real_t volume = mesh.get_volume(cell_id_old);

    real_t s = (particle.t - t_old) * weight_old / volume * particle.intensity;

    density[cell_id_old] += s;

    momentum_x[cell_id_old] += particle_old.geometry.velocity.x * s;
    momentum_y[cell_id_old] += particle_old.geometry.velocity.y * s;
    momentum_z[cell_id_old] += particle_old.geometry.velocity.z * s;

    energy[cell_id_old] +=
        s * (0.5 * particle_old.mass * particle_old.geometry.velocity.norm2()) *
        constants::m_p / constants::eV_to_Joules * (1/W);

    energy_x[cell_id_old] +=
        s * (0.5 * particle_old.mass * particle_old.geometry.velocity.x * particle_old.geometry.velocity.x) *
        constants::m_p / constants::eV_to_Joules * (1/W);
}

void ArcmcEstimator::add_final(const Particle& particle, const Mesh& mesh, real_t W)
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
        (0.5 * particle.mass * particle.geometry.velocity.norm2()) *
        particle.weight / volume * particle.intensity * (1/W);

    energy_x_final[cell_id] +=
        (0.5 * particle.mass * particle.geometry.velocity.x * particle.geometry.velocity.x) *
        particle.weight / volume * particle.intensity * (1/W);
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