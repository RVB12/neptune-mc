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

    // Only for voxel meshes, also check if the particle is likely to cross
    // the boundary
    const VoxelMesh* voxel_mesh = dynamic_cast<const VoxelMesh*>(&mesh);
    if (voxel_mesh != nullptr) {
        Vec3 temp_diff(D, D, D);
        Vec3 temp_min = p.geometry.position + drift_term - 2 * temp_diff;
        Vec3 temp_max = p.geometry.position + drift_term + 2 * temp_diff;

        if (temp_min.x < 0 or temp_max.x > voxel_mesh->Lx or temp_min.y < 0 or
            temp_max.y > voxel_mesh->Ly or temp_min.z < 0 or
            temp_max.z > voxel_mesh->Lz) {
            return false;
        }
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

// Generate a sample from N(0, 1) truncated to [a, ∞) using Robert's algorithm
// Note, this will get more and more expensive for larger a, it may
// need to be replaced with a fancier algorithm at some point
real_t truncated_normal_right_tail(real_t a, std::mt19937& rng)
{
    if (a <= 0.0) {
        // For small a, use rejection from standard normal
        std::normal_distribution<> normal(0.0, 1.0);
        real_t x;
        do {
            x = normal(rng);
        } while (x < a);
        return x;
    }

    // For a > 0, use exponential proposal
    real_t alpha = (a + std::sqrt(a * a + 4.0)) / 2.0;
    std::exponential_distribution<> exponential(alpha);
    std::uniform_real_distribution<> uniform(0.0, 1.0);

    while (true) {
        real_t z = exponential(rng) + a;
        real_t u = uniform(rng);
        if (u <= std::exp(-0.5 * std::pow(z - alpha, 2))) {
            return z;
        }
    }
}

real_t normal_cdf(real_t x, real_t mean, real_t variance)
{
    return 0.5 * (1.0 + std::erf((x - mean) / std::sqrt(2.0 * variance)));
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

    real_t D0 = bgk_stddev * bgk_stddev / (R_cx * R_cx_theta);
    D0 *= (std::exp(-R_cx_theta) - 1 + R_cx_theta);
    if (D0 < 0) {
        D0 = 0;
    }
    real_t D = std::sqrt(2 * D0 * dt);

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

                Vec3 normal = mesh.get_surface_normal(p.geometry.cell_id,
                                                      crossing.surface_id);

                // Make sure the normal points outwards
                if (dot(normal, p.geometry.velocity) < 0) {
                    // The normal points inwards, we need to flip it
                    normal = -normal;
                }

                real_t drift_normal = dot(u, normal);
                if (std::abs(drift_normal) < 1.0e-6) {
                    // There is no (significant) drift in the direction
                    // of the normal, simplified case

                    // Velocity can just be reflected
                    p.geometry.velocity =
                        p.geometry.velocity -
                        2 * dot(p.geometry.velocity, normal) * normal;

                    AbsorbingBC* absorbing_bc = dynamic_cast<AbsorbingBC*>(
                        boundary_conditions[crossing.boundary_condition_id]
                            .get());
                    if (absorbing_bc != nullptr) {
                        p.weight *= -1;
                    }
                    // TODO: handle partial absorption
                }
                else {
                    // There is a drift, we can't just reflect the velocity,
                    // we need to resample the new position from the tail of
                    // the mirrored transition function
                    // potentially even sampling from an infinite amount
                    // of reflected solutions

                    AbsorbingBC* absorbing_bc = dynamic_cast<AbsorbingBC*>(
                        boundary_conditions[crossing.boundary_condition_id]
                            .get());
                    if (absorbing_bc != nullptr) {
                        // Absorbing is the "easy" case,
                        // we only have one truncated gaussian to sample from

                        Vec3 mirror_pos =
                            old_pos -
                            2 * dot(old_pos - crossing.position, normal) *
                                normal;
                        Vec3 new_center = mirror_pos + u * dt;

                        // Get the two vectors that are orthogonal to the normal
                        // and to each other using Gram-Schmidt
                        Vec3 v2(r_normal(p.rng), r_normal(p.rng),
                                r_normal(p.rng));
                        v2 -= dot(v2, normal) * normal;
                        v2 /= v2.norm();
                        Vec3 v3 = normal.cross(v2);

                        real_t n2 = r_normal(p.rng);
                        real_t n3 = r_normal(p.rng);

                        // Sample n1 from the normal distribution constrained to
                        // n1
                        real_t dist_to_crossing =
                            -dot(crossing.position - new_center, normal);
                        real_t n1 = -truncated_normal_right_tail(
                            dist_to_crossing / D, p.rng);

                        Vec3 new_pos2 = new_center + n1 * D * normal +
                                        n2 * D * v2 + n3 * D * v3;

                        // This is an extra check that should not really be
                        // necessary, check if the new position is actually
                        // inside the domain
                        if (dot(new_pos2 - crossing.position, normal) >= 0) {
                            // We are inside the cell, we can break
                            throw std::runtime_error(
                                "New position is not inside the domain");
                        }

                        // Determine a new assumed velocity
                        // so that the particle will move to the new position
                        p.geometry.velocity =
                            (1 / time_left) * (new_pos2 - p.geometry.position);

                        //
                        // Now we need to update the weight:
                        //

                        Vec3 old_center = old_pos + u * dt;
                        // Calculate area from the original distribution
                        // that is outside of the domain
                        real_t distance_left =
                            dot(crossing.position - old_center, normal);
                        real_t area_left =
                            1 - normal_cdf(distance_left, 0, 2 * D0 * dt);
                        // Calculate area from the new distribution
                        // that is inside the domain
                        real_t distance_right =
                            dot(crossing.position - new_center, normal);
                        real_t area_right =
                            normal_cdf(distance_right, 0, 2 * D0 * dt);
                        real_t A =
                            -std::exp(drift_normal / D0 *
                                      dot(crossing.position - old_pos, normal));

                        p.weight *= A * area_right / area_left;

                        old_pos = mirror_pos;
                    }
                    else {
                        // reflecting boundary condition
                        // the "hard" case
                        // we either need to sample from a normal
                        // centered around the mirror position
                        // or we need to sample from a position
                        // along a line of positions.

                        Vec3 mirror_pos =
                            old_pos -
                            2 * dot(old_pos - crossing.position, normal) *
                                normal;
                        Vec3 old_center = old_pos + u * dt;
                        Vec3 new_center = mirror_pos + u * dt;

                        // First, we determine whether we need to sample
                        // from the point source or from the line
                        // for that, we need their total area

                        // The area of the point source is easy,
                        // it is just a normal distribution
                        real_t d1 = dot(crossing.position - new_center, normal);
                        real_t area_point = normal_cdf(d1, 0, 2 * D0 * dt);

                        real_t alpha = dot(u, normal);
                        real_t beta = D0;
                        real_t A =
                            std::exp(drift_normal / D0 *
                                     dot(crossing.position - old_pos, normal));

                        area_point *= A;

                        // The area of the line is more complicated,
                        // we use the property that for pure reflection,
                        // the areas should sum to 1
                        real_t d2 = dot(crossing.position - old_center, normal);
                        real_t area_original =
                            1 - normal_cdf(d2, 0, 2 * D0 * dt) - area_point;

                        real_t u1 = r_unif(p.rng);
                        if (u1 < area_point / (area_point + area_original)) {
                            // Sample from the point source
                        }
                        else {
                            while (true) {
                                // Sample from the line
                                std::exponential_distribution<real_t>
                                    exponential(alpha / beta);

                                real_t offset = exponential(p.rng);

                                // // TO BE CHECKED: do i need to correct for
                                // the fact
                                // // that the area of the truncated normal is
                                // not 1?

                                real_t area_new_point =
                                    normal_cdf(d1 + offset, 0, 2 * D0 * dt);
                                real_t u2 = r_unif(p.rng);
                                if (u2 < area_new_point / area_point) {
                                    new_center =
                                        new_center +
                                        offset *
                                            normal;  // Move along the normal
                                    break;
                                }
                            }
                        }

                        // Now we can sample from the shifted normal

                        // Get the two vectors that are orthogonal to the normal
                        // and to each other using Gram-Schmidt
                        Vec3 v2(r_normal(p.rng), r_normal(p.rng),
                                r_normal(p.rng));
                        v2 -= dot(v2, normal) * normal;
                        v2 /= v2.norm();
                        Vec3 v3 = normal.cross(v2);

                        real_t n2 = r_normal(p.rng);
                        real_t n3 = r_normal(p.rng);

                        // Sample n1 from the normal distribution constrained to
                        // n1
                        real_t dist_to_crossing =
                            -dot(crossing.position - new_center, normal);
                        real_t n1 = -truncated_normal_right_tail(
                            dist_to_crossing / D, p.rng);

                        Vec3 new_pos2 = new_center + n1 * D * normal +
                                        n2 * D * v2 + n3 * D * v3;

                        // This is an extra check that should not really be
                        // necessary, check if the new position is actually
                        // inside the domain
                        if (dot(new_pos2 - crossing.position, normal) >= 0) {
                            // We are inside the cell, we can break
                            throw std::runtime_error(
                                "New position is not inside the domain");
                        }

                        // Determine a new assumed velocity
                        // so that the particle will move to the new position
                        p.geometry.velocity =
                            (1 / time_left) * (new_pos2 - p.geometry.position);

                        old_pos = mirror_pos;
                    }
                }

                // TODO: handle partial absorption
            }
        }
        else {
            // We will not hit a cell boundary before the time is up
            // Move the particle to the new position and draw a new velocity
            p.geometry.position =
                p.geometry.position + time_left * p.geometry.velocity;
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