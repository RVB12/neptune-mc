#include <omp.h>

#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/linesource.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/physics/update_particle.hpp"
#include "neptune/simulator/simulator.hpp"
namespace neptune {

/**
 * @brief Simulate a particle using a diffusive fluid approximation
 *
 * We assume there is only one bgk process
 *
 * This is the old version of the algorithm by Bert
 *
 *
 * @return true if the step was successful
 *         false if the particle left the mesh, in that case, we must
 *         switch to a kinetic simulation
 */
bool kdmc_fluid_step_old(
    const Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    Particle& p,
    real_t dt)
{
    if (dt <= 0) {
        return true;
    }
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    real_t R_cx = bgk_processes[0]->rate[cell_id];
    real_t R_cx_theta = R_cx * dt;

    Vec3 u = Vec3(bgk_processes[0]->bgk_mean_x[cell_id],
                  bgk_processes[0]->bgk_mean_y[cell_id],
                  bgk_processes[0]->bgk_mean_z[cell_id]);

    Vec3 drift_term = u * dt;

    real_t bgk_stddev = bgk_processes[0]->bgk_stddev[cell_id];

    real_t D = bgk_stddev * bgk_stddev / (R_cx * R_cx_theta);
    D *= (std::exp(-R_cx_theta) - 1 + R_cx_theta);
    if (D < 0) {
        D = 0;
    }
    D = std::sqrt(2 * D * dt);
    // TODO, add velocity induced diffusion and hetergeneity correction

    Vec3 diffusion_term;
    diffusion_term.x = D * r_normal(p.rng);
    diffusion_term.y = D * r_normal(p.rng);
    diffusion_term.z = D * r_normal(p.rng);

    Vec3 new_pos = p.geometry.position + drift_term + diffusion_term;

    // Check if the new position is inside the mesh
    index_t new_cell_id = mesh.get_cell_index(new_pos);
    if (new_cell_id == -1) {
        return false;
    }

    // Update the particle
    p.geometry.position = new_pos;
    p.geometry.cell_id = new_cell_id;
    p.geometry.boundary_id = -1;
    p.t += dt;
    if (ion_processes.size() > 0) {
        p.weight *= std::exp(-dt * ion_processes[0]->rate[cell_id]);
    }
    return true;
}

/**
 * @brief Simulate a particle using a diffusive fluid approximation
 *
 * We assume there is only one bgk process
 */
void kdmc_fluid_step(
    const Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    Particle& p,
    real_t dt)
{
    if (dt <= 0) {
        return;
    }
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    real_t R_cx = bgk_processes[0]->rate[cell_id];
    real_t R_cx_theta = R_cx * dt;

    Vec3 u = Vec3(bgk_processes[0]->bgk_mean_x[cell_id],
                  bgk_processes[0]->bgk_mean_y[cell_id],
                  bgk_processes[0]->bgk_mean_z[cell_id]);

    Vec3 drift_term = u * dt;

    real_t bgk_stddev = bgk_processes[0]->bgk_stddev[cell_id];

    real_t D = bgk_stddev * bgk_stddev / (R_cx * R_cx_theta);
    D *= (std::exp(-R_cx_theta) - 1 + R_cx_theta);
    if (D < 0) {
        D = 0;
    }
    D = std::sqrt(2 * D * dt);

    Vec3 diffusion_term;
    diffusion_term.x = D * r_normal(p.rng);
    diffusion_term.y = D * r_normal(p.rng);
    diffusion_term.z = D * r_normal(p.rng);

    Vec3 old_pos = p.geometry.position;
    Vec3 new_pos = old_pos + drift_term + diffusion_term;

    // Pretend like the particle is moving in a straight line from the old the
    // the new pos
    Vec3 assumed_velocity = (1 / dt) * (new_pos - p.geometry.position);
    p.geometry.velocity = assumed_velocity;

    real_t time_left = dt;
    while (true) {
        CellCrossing crossing = mesh.get_next_crossing(p.geometry);

        if (crossing.t < time_left) {
            // We will hit a boundary before the time is up
            // Move the particle to the boundary
            p.t += crossing.t;
            time_left -= crossing.t;
            p.geometry.position = crossing.position;

            if (crossing.next_cell_id != -1) {
                // Regular grid crossing, move the particle to the next
                // cell
                p.geometry.boundary_id = crossing.next_surface_id;
                p.geometry.cell_id = crossing.next_cell_id;
            }
            else {
                // Boundary crossing
                p.geometry.boundary_id = crossing.surface_id;

                Vec3 normal =
                    mesh.get_surface_normal(cell_id, crossing.surface_id);

                real_t A = 0;
                SpecularBC* specular_bc = dynamic_cast<SpecularBC*>(
                    boundary_conditions[crossing.surface_id].get());
                if (specular_bc != nullptr) {
                    A = 1;
                }
                else {
                    AbsorbingBC* absorbing_bc = dynamic_cast<AbsorbingBC*>(
                        boundary_conditions[crossing.surface_id].get());
                    if (absorbing_bc != nullptr) {
                        // Note that for an absorbing boundary, the weight is
                        // not set to zero, but is multiplied by -1
                        A = -1;
                    }
                    else {
                        throw std::runtime_error(
                            "Using a boundary condition that is not supported "
                            "by KD");
                    }
                }
                p.weight *= A;

                // Reflect the particle
                if (u.x == 0 && u.y == 0 && u.z == 0) {
                    // No drift, we can just reflect the velocity
                    p.geometry.velocity =
                        p.geometry.velocity -
                        2 * (p.geometry.velocity * normal) * normal;
                }
                else {
                    // There is a drift, we can't just reflect the velocity,
                    // we need to resample the new position from the tail of
                    // the mirrored transition function

                    // TODO: make this more general, this only works for the furthest boundary
                    const VoxelMesh* mesh_ptr =
                        dynamic_cast<const VoxelMesh*>(&mesh);

                    real_t Lx = mesh_ptr->Lx;

                    Vec3 mirror_pos = old_pos;
                    mirror_pos.x = 2 * Lx - mirror_pos.x;

                    Vec3 new_center = mirror_pos + u * dt;

                    Vec3 new_pos2;
                    while(true) {
                        // Sample a new position from the mirror distribution
                        Vec3 diffusion_term2 = Vec3(
                            D * r_normal(p.rng), D * r_normal(p.rng),
                            D * r_normal(p.rng));
                        new_pos2 = new_center + diffusion_term2;

                        if(new_pos2.x < Lx){
                            break;
                        }
                    }

                    // Determine a new assumed velocity
                    p.geometry.velocity = (1/time_left) * (new_pos2 - p.geometry.position);

                }
            }
        }
        else {
            // We will not hit a cell boundary before the time is up
            // Move the particle to the new position and draw a new velocity
            p.geometry.position = new_pos;
            p.geometry.boundary_id = -1;
            p.t += time_left;
            time_left = 0;

            index_t cell_id = p.geometry.cell_id;
            // Draw a new velocity
            std::normal_distribution<real_t> dist_vx(
                bgk_processes[0]->bgk_mean_x[cell_id], bgk_stddev);
            std::normal_distribution<real_t> dist_vy(
                bgk_processes[0]->bgk_mean_y[cell_id], bgk_stddev);
            std::normal_distribution<real_t> dist_vz(
                bgk_processes[0]->bgk_mean_z[cell_id], bgk_stddev);
            p.geometry.velocity =
                Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));
            break;
        }
    }

    if (ion_processes.size() > 0) {
        p.weight *= std::exp(-dt * ion_processes[0]->rate[cell_id]);
    }
}

void simulate_kdmc(
    Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    vector<std::unique_ptr<Source>>& sources,
    Estimator& estimator,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats,
    real_t dt,
    bool use_new_version)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // History based simulation (because its easy)
#pragma omp parallel
    {
        index_t num_crossings_local = 0;
        index_t num_collisions_local = 0;

        Estimator estimator_local(mesh.n_cells());
        for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp for
            for (int i_particle = 0;
                 i_particle < sources[i_source]->num_particles();
                 i_particle++) {
                // Generate particle
                Particle particle =
                    sources[i_source]->generate_particle(i_particle);
                particle.geometry.cell_id =
                    mesh.get_cell_index(particle.geometry.position);

                particle.intensity = sources[i_source]->intensity() /
                                     sources[i_source]->num_particles();

                if (particle.geometry.cell_id == -1) {
                    throw std::runtime_error(
                        "Particle generated outside of mesh");
                }

                for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                    Particle particle_old = particle;

                    // Simulate particle kinetically until a charge exchange
                    // collision occurs

                    real_t t_max = particle.t + dt;
                    if (params.t_final > 0) {
                        t_max = std::min(t_max, params.t_final);
                    }

                    while (true) {
                        particle_old = particle;
                        EventType event = update_particle(
                            mesh, bgk_processes, ion_processes,
                            boundary_conditions, particle, t_max);

                        estimator_local.add_track_length(particle, particle_old,
                                                         mesh);
                        if (particle.weight < 1.0e-4) {
                            break;
                        }
                        if (event == EventType::FINAL_TIME) {
                            break;
                        }
                        else if (event == EventType::CHARGE_EXCHANGE) {
                            num_collisions_local++;
                            break;
                        }
                        else if (event == EventType::CELL_CROSSING) {
                            num_crossings_local++;
                        }
                    }

                    // Now simulate the particle fluidly for the rest of the
                    // time step

                    if (particle.t < t_max) {
                        // Now determine the final free flight time

                        real_t final_free_flight_time;

                        if (use_new_version) {
                            // number of collisions is a poisson distribution
                            std::poisson_distribution<index_t> poisson(
                                bgk_processes[0]
                                    ->rate[particle.geometry.cell_id] *
                                (t_max - particle.t));
                            index_t n_collisions = poisson(particle.rng);
                            // now that we have the number of collisions, the
                            // time of the last event is a beta distribution C++
                            // does not have a beta distribution, so we use the
                            // two gamma distribution trick
                            std::gamma_distribution<real_t> gamma1(n_collisions,
                                                                   1);
                            std::gamma_distribution<real_t> gamma2(1, 1);
                            real_t X = gamma1(particle.rng);
                            real_t Y = gamma2(particle.rng);
                            final_free_flight_time =
                                (1 - X / (X + Y)) * (t_max - particle.t);
                        }
                        else {
                            real_t R_cx = bgk_processes[0]
                                              ->rate[particle.geometry.cell_id];
                            real_t eps = r_unif(particle.rng);
                            final_free_flight_time = -std::log(eps) / R_cx;
                        }

                        if (final_free_flight_time > t_max - particle.t) {
                            final_free_flight_time = t_max - particle.t;
                        }
                        else {
                            if (use_new_version) {
                                kdmc_fluid_step(mesh, bgk_processes,
                                                ion_processes,
                                                boundary_conditions, particle,
                                                t_max - particle.t -
                                                    final_free_flight_time);

                                // TODO: add fluid estimator here
                            }
                            else {
                                bool step_succesfull = kdmc_fluid_step_old(
                                    mesh, bgk_processes, ion_processes,
                                    boundary_conditions, particle,
                                    t_max - particle.t -
                                        final_free_flight_time);

                                if (!step_succesfull) {
                                    // Particle left the mesh, switch to kinetic
                                    // simulation
                                    final_free_flight_time = t_max - particle.t;
                                    while (true) {
                                        particle_old = particle;
                                        EventType event = update_particle(
                                            mesh, bgk_processes, ion_processes,
                                            boundary_conditions, particle,
                                            t_max);

                                        estimator_local.add_track_length(
                                            particle, particle_old, mesh);
                                        if (particle.weight < 1.0e-4) {
                                            break;
                                        }
                                        if (event == EventType::FINAL_TIME) {
                                            break;
                                        }
                                        else if (event ==
                                                 EventType::CHARGE_EXCHANGE) {
                                            num_collisions_local++;
                                        }
                                        else if (event ==
                                                 EventType::CELL_CROSSING) {
                                            num_crossings_local++;
                                        }
                                    }
                                }
                                else {
                                    // TODO: add fluid estimator here
                                }
                            }
                        }

                        // Final free flight time
                        while (particle.t < t_max) {
                            particle_old = particle;
                            // Create empty collision lists to simulate no
                            // collisions
                            vector<std::unique_ptr<BGKScatteringProcess>>
                                empty_bgk_processes;
                            vector<std::unique_ptr<IonisationProcess>>
                                empty_ion_processes;
                            EventType event = update_particle(
                                mesh, empty_bgk_processes, empty_ion_processes,
                                boundary_conditions, particle, t_max);

                            estimator_local.add_track_length(
                                particle, particle_old, mesh);
                            if (particle.weight < 1.0e-4) {
                                break;
                            }
                            if (event == EventType::FINAL_TIME) {
                                break;
                            }
                            else if (event == EventType::CHARGE_EXCHANGE) {
                                num_collisions_local++;
                            }
                            else if (event == EventType::CELL_CROSSING) {
                                num_crossings_local++;
                            }
                        }
                    }

                    if (std::abs(particle.weight) < 1.0e-4) {
                        break;
                    }

                    // Tally
                    if (particle.t >= params.t_final) {
                        estimator_local.add_final(particle, mesh);
                        break;
                    }
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

}  // namespace neptune