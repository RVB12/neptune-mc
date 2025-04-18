#include <omp.h>

#include "eigrene/geometry/voxelmesh.hpp"
#include "eigrene/physics/linesource.hpp"
#include "eigrene/physics/physicalconstants.hpp"
#include "eigrene/physics/update_particle.hpp"
#include "eigrene/simulator/simulator.hpp"
namespace eigrene {

void update_particle_kdmc2(const Mesh& mesh,
                           const PlasmaBackground& background,
                           const vector<BoundaryCondition>& boundary_conditions,
                           Particle& p,
                           index_t& num_crossings,
                           index_t& num_collisions,
                           real_t t_final,
                           real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    // Limit the time step to the final time
    real_t dt2 = std::min(dt, t_final - p.t);

    //
    // Step 1: Kinetic simulation until a collision occurs
    //

    real_t t_kinetic = 0;
    while (true) {
        // random number to determine time till collision
        real_t eps = r_unif(p.rng);
        // time till collision
        real_t tc =
            -std::log(eps) /
            (constants::calculate_K_cx(background.plasma_temp[cell_id]) *
             background.plasma_dens[cell_id]);

        CellCrossing crossing = mesh.get_next_crossing(p.geometry);

        // Check which one is earlier, charge exchange, cell crossing or final
        // time
        real_t time_left = dt2 - t_kinetic;
        if (time_left < tc && time_left < crossing.t) {
            // We will not hit the boundary before the time is up
            // This is the final iteration
            p.geometry.position += time_left * p.geometry.velocity;
            p.geometry.boundary_id = -1;
            t_kinetic = dt2;
            p.weight =
                p.weight * std::exp(-time_left *
                                    constants::calculate_K_i(
                                        background.plasma_temp[cell_id]) *
                                    background.plasma_dens[cell_id]);
            break;
        }
        if (tc < crossing.t) {
            num_collisions++;
            // Move the particle until the collision
            p.geometry.position += tc * p.geometry.velocity;
            // Generate new velocity from plasma background (BGK
            // operator)
            real_t stddev = constants::calculate_BGK_stddev(
                background.plasma_temp[cell_id]);
            std::normal_distribution<real_t> dist_vx(
                background.plasma_vel_x[cell_id], stddev);
            std::normal_distribution<real_t> dist_vy(
                background.plasma_vel_y[cell_id], stddev);
            std::normal_distribution<real_t> dist_vz(
                background.plasma_vel_z[cell_id], stddev);
            p.geometry.velocity =
                Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

            // Reset boundary id (we are no longer at a boundary)
            p.geometry.boundary_id = -1;
            t_kinetic += tc;
            // We've had a collision, time for the diffusive step
            break;
        }
        else {
            num_crossings++;
            t_kinetic += crossing.t;
            if (crossing.next_cell_id != -1) {
                // Regular grid crossing, move the particle to the next
                // cell
                p.geometry.position = crossing.position;
                p.geometry.cell_id = crossing.next_cell_id;
                p.geometry.boundary_id = crossing.next_surface_id;
            }
            else {
                // Boundary crossing
                p.geometry.position = crossing.position;
                p.geometry.boundary_id = crossing.surface_id;
                BoundaryCondition bc =
                    boundary_conditions[crossing.boundary_condition_id];

                if (bc == BoundaryCondition::REFLECTIVE) {
                    // Reflect the particle
                    Vec3 normal =
                        mesh.get_surface_normal(cell_id, crossing.surface_id);
                    p.geometry.velocity =
                        p.geometry.velocity -
                        2 * dot(p.geometry.velocity, normal) * normal;
                }
                else if (bc == BoundaryCondition::ABSORBING) {
                    // Absorb the particle
                    p.weight = 0;
                    break;
                }
                else {
                    throw std::runtime_error(
                        "Using a boundary condition that is not (yet) "
                        "implemented");
                }
            }
        }
    }

    if (t_kinetic >= dt2 or std::abs(p.weight) <= 0) {
        // We've reached the final time, no need to do the diffusive step
        p.t += t_kinetic;
        return;
    }

    //
    // Step 2: Diffusive simulation for the rest of the time step
    //

    // The particle will not actually move diffusively for the
    // entire diffusive time. We will also take the time between
    // the final collision in the diffusive step and the end of the
    // time step into account.
    //
    // In a homogeneous medium, we can sample this as
    // 1. Draw the number of collisions from a Poisson distribution
    // 2. Draw the time of the final collision from a beta distribution
    // In a heterogeneous medium, this gets more complicated, but
    // we will just use the time from the homogeneous case.
    real_t t_diffusive = dt2 - t_kinetic;

    real_t Rcx = constants::calculate_K_cx(background.plasma_temp[cell_id]) *
                 background.plasma_dens[cell_id];

    // Commented out for now until final kinetic time is fixed
    // std::poisson_distribution<int> poisson_dist(Rcx * t_diffusive);
    // int num_collisions = poisson_dist(p.rng);
    // {
    //     // Draw sample from beta distribution
    //     // C++ doesn't have a built-in beta distribution, so we will
    //     // use two gamma distributions
    //     std::gamma_distribution<real_t> gamma_dist1(num_collisions, 1);
    //     std::gamma_distribution<real_t> gamma_dist2(1, 1);

    //     real_t gamma1 = gamma_dist1(p.rng);
    //     real_t gamma2 = gamma_dist2(p.rng);
    //     t_diffusive = gamma1 / (gamma1 + gamma2) * t_diffusive;
    // }

    //
    // During a diffusive step, the particle behaves fluidly according to the
    // SDE: dX = u(X) dt + D(X) dW For now, we just take one large timestep
    // (this is exact for homogeneous grids) but in the future we should look
    // into higher order schemes and/or smaller timesteps to sample from this
    // SDE.
    //

    // The diffusion coefficient is not just the expected diffusion of the fluid
    // we add an exponential term to that for long timesteps, it converges to
    // the diffusive limit, but converges to kinetic behavior for small
    // timesteps
    real_t Rcx_theta = Rcx * t_diffusive;
    real_t D = background.plasma_temp[cell_id] * constants::eV_to_Joules /
               constants::m_p / (Rcx * Rcx_theta);
    D *= (std::exp(-Rcx_theta) - 1 + Rcx_theta);
    if (D < 0) {
        D = 0;
    }

    // The drift is equal to the velocity of the plasma background
    Vec3 u =
        Vec3(background.plasma_vel_x[cell_id], background.plasma_vel_y[cell_id],
             background.plasma_vel_z[cell_id]);

    // Sample brownian step
    Vec3 W = Vec3(r_normal(p.rng), r_normal(p.rng), r_normal(p.rng)) *
             std::sqrt(t_diffusive);

    // Compute the new position
    Vec3 new_position =
        p.geometry.position + u * t_diffusive + std::sqrt(2 * D) * W;
    index_t new_cell_id = mesh.get_cell_index(new_position);

    // Now handle boundary conditions
    // Assume the boundary conditions are the same everywhere
    if (boundary_conditions[0] == BoundaryCondition::REFLECTIVE) {
        if (new_cell_id == -1) {
            VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);

            // Reflect the particle back inside the mesh
            // This may have to be done multiple times
            while (true) {
                if (new_position.x < 0) {
                    new_position.x = -new_position.x;
                }
                else if (new_position.x > vmesh.Lx) {
                    new_position.x = 2 * vmesh.Lx - new_position.x;
                }
                if (new_position.y < 0) {
                    new_position.y = -new_position.y;
                }
                else if (new_position.y > vmesh.Ly) {
                    new_position.y = 2 * vmesh.Ly - new_position.y;
                }
                if (new_position.z < 0) {
                    new_position.z = -new_position.z;
                }
                else if (new_position.z > vmesh.Lz) {
                    new_position.z = 2 * vmesh.Lz - new_position.z;
                }

                new_cell_id = mesh.get_cell_index(new_position);
                if (new_cell_id != -1) {
                    break;
                }
            }
        }
    }
    else {
        if (new_cell_id == -1) {
            VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);

            // Reflect the particle back inside the mesh
            // This may have to be done multiple times
            while (true) {
                if (new_position.x < 0) {
                    new_position.x = -new_position.x;
                }
                else if (new_position.x > vmesh.Lx) {
                    new_position.x = 2 * vmesh.Lx - new_position.x;
                }
                if (new_position.y < 0) {
                    new_position.y = -new_position.y;
                }
                else if (new_position.y > vmesh.Ly) {
                    new_position.y = 2 * vmesh.Ly - new_position.y;
                }
                if (new_position.z < 0) {
                    new_position.z = -new_position.z;
                }
                else if (new_position.z > vmesh.Lz) {
                    new_position.z = 2 * vmesh.Lz - new_position.z;
                }

                new_cell_id = mesh.get_cell_index(new_position);
                if (new_cell_id != -1) {
                    break;
                }
            }

            // make the weight negative, not just zero
            p.weight = -p.weight;
        }
        // Absorbing boundary condition
        // if (new_cell_id == -1) {
        //     p.weight = 0;
        //     p.t += dt2;
        //     return;
        // } else {
        //     // Even if the new position is inside the mesh
        //     // it could have crossed the boundary and gone back
        //     // in the diffusive time step. We need to correct for that
        //     // Using the method of images, the new pdf would be the
        //     unconstrained
        //     // pdf - the reflected pdf. We have sampled from the original pdf
        //     // so we need to subtract all the reflected pdfs and do
        //     accept-reject

        //     // Note: these equations still need to be fixed for drift cases
        //     // and also for other boundaries

        //     real_t x_center = p.geometry.position.x;

        //     real_t x_new = new_position.x;

        //     real_t p_original =
        //         std::exp(-(x_new - x_center)*(x_new - x_center)/(2 * D *
        //         t_diffusive)) / std::sqrt(4 * M_PI * D * t_diffusive);

        //     VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);

        //     real_t x_image = 2 * vmesh.Lx - x_center;

        //     real_t p_reflected =
        //         std::exp(-(x_new - x_image)*(x_new - x_image)/(2 * D *
        //         t_diffusive)) / std::sqrt(4 * M_PI * D * t_diffusive);
        //     real_t p_total = p_original - p_reflected;

        //     // std::cout << "p_original: " << p_original << std::endl;
        //     // std::cout << "p_reflected: " << p_reflected << std::endl;
        //     // std::cout << "p_total: " << p_total << std::endl;

        //     // real_t u = r_unif(p.rng);
        //     // if( u * p_original > p_total) {
        //     //     // We have to reject the sample
        //     //     p.weight = 0;
        //     //     p.t += dt2;
        //     //     return;
        //     // }

        // }
    }

    p.geometry.position = new_position;
    p.geometry.cell_id = new_cell_id;

    //
    // Step 3: Final free flight time
    //

    real_t t_kinetic2 = dt2 - t_kinetic - t_diffusive;

    // First, sample a new velocity

    real_t stddev =
        constants::calculate_BGK_stddev(background.plasma_temp[new_cell_id]);
    std::normal_distribution<real_t> dist_vx(
        background.plasma_vel_x[new_cell_id], stddev);
    std::normal_distribution<real_t> dist_vy(
        background.plasma_vel_y[new_cell_id], stddev);
    std::normal_distribution<real_t> dist_vz(
        background.plasma_vel_z[new_cell_id], stddev);

    p.geometry.velocity = Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

    // // Now simulate the particle kinetically for the rest of the time step
    // // However, we don't need to check for charge exchange as we have already
    // // determined that no charge exchange will occur.

    p.t += dt2;
}

void simulate_kdmc2(Mesh& mesh,
                    PlasmaBackground& background,
                    vector<std::unique_ptr<Source>>& sources,
                    Estimator& estimator,
                    vector<BoundaryCondition>& boundary_conditions,
                    SimulationParams params,
                    ExecutionStatistics& stats,
                    real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // History based simulation (because its easy)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp parallel
        {
            index_t num_crossings_local = 0;
            index_t num_collisions_local = 0;

            Estimator estimator_local(mesh.n_cells());
#pragma omp for
            for (int i_particle = 0;
                 i_particle < sources[i_source]->num_particles();
                 i_particle++) {
                real_t intensity_multiplier =
                    sources[i_source]->intensity() /
                    sources[i_source]->num_particles();

                // Generate particle
                Particle particle =
                    sources[i_source]->generate_particle(i_particle);
                particle.geometry.cell_id =
                    mesh.get_cell_index(particle.geometry.position);

                if (particle.geometry.cell_id == -1) {
                    throw std::runtime_error(
                        "Particle generated outside of mesh");
                }

                for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                    Particle particle_old = particle;
                    // Simulate the particle up to the next event
                    update_particle_kdmc2(mesh, background, boundary_conditions,
                                          particle, num_crossings_local,
                                          num_collisions_local, params.t_final,
                                          dt);

                    // TODO: add fluid estimator
                    estimator_local.add_track_length(particle, particle_old,
                                                     mesh);

                    if (std::abs(particle.weight) < 1.0e-4) {
                        break;
                    }

                    // Tally
                    if (particle.t >= params.t_final) {
                        real_t v = mesh.get_volume(particle.geometry.cell_id);
                        estimator_local.add_final(particle, mesh);
                        break;
                    }
                }
            }

#pragma omp critical
            {
                // Add the local results to the global results
                estimator.combine(estimator_local);
            }
#pragma omp atomic
            stats.num_crossings += num_crossings_local;
#pragma omp atomic
            stats.num_collisions += num_collisions_local;
        }
    }
}

}  // namespace eigrene