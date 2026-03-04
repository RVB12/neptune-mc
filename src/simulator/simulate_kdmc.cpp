#include <omp.h>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/physics/update_particle.hpp"
#include "neptune/simulator/simulator.hpp"
namespace neptune {

real_t normal_cdf(real_t x, real_t mean, real_t variance)
{
    return 0.5 * (1.0 + std::erf((x - mean) / std::sqrt(2.0 * variance)));
}

Vec3 calculate_drift(
    const Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    real_t dt,
    index_t cell_id)
{
    Vec3 u = Vec3(bgk_processes[0]->bgk_mean_x[cell_id],
                  bgk_processes[0]->bgk_mean_y[cell_id],
                  bgk_processes[0]->bgk_mean_z[cell_id]);

    real_t R_cx = bgk_processes[0]->rate[cell_id];
    real_t R_cx_theta = R_cx * dt;

    real_t bgk_stddev = bgk_processes[0]->bgk_stddev[cell_id];

    // Calculate derivative of the diffusion term and add it to the drift term
    Vec3 diffusion_induced_drift(0, 0, 0);
    VoxelMesh* voxel_mesh = dynamic_cast<VoxelMesh*>(const_cast<Mesh*>(&mesh));
    if (voxel_mesh != nullptr) {
        Index3D ij = voxel_mesh->linear_to_tuple(cell_id);

        if (ij.i > 0 and ij.i < voxel_mesh->nx - 1) {
            diffusion_induced_drift.x =
                (bgk_processes[0]->rate[voxel_mesh->tuple_to_linear(
                     Index3D(ij.i + 1, ij.j, ij.k))] -
                 bgk_processes[0]->rate[voxel_mesh->tuple_to_linear(
                     Index3D(ij.i - 1, ij.j, ij.k))]) /
                (2 * voxel_mesh->dx);
        }

        if (ij.j > 0 and ij.j < voxel_mesh->ny - 1) {
            diffusion_induced_drift.y =
                (bgk_processes[0]->rate[voxel_mesh->tuple_to_linear(
                     Index3D(ij.i, ij.j + 1, ij.k))] -
                 bgk_processes[0]->rate[voxel_mesh->tuple_to_linear(
                     Index3D(ij.i, ij.j - 1, ij.k))]) /
                (2 * voxel_mesh->dy);
        }

        if (ij.k > 0 and ij.k < voxel_mesh->nz - 1) {
            diffusion_induced_drift.z =
                (bgk_processes[0]->rate[voxel_mesh->tuple_to_linear(
                     Index3D(ij.i, ij.j, ij.k + 1))] -
                 bgk_processes[0]->rate[voxel_mesh->tuple_to_linear(
                     Index3D(ij.i, ij.j, ij.k - 1))]) /
                (2 * voxel_mesh->dz);
        }

        real_t dDdsigma = (-2 * bgk_stddev * bgk_stddev) /
                          (R_cx * R_cx * R_cx_theta) *
                          (std::exp(-R_cx_theta) - 1 + R_cx_theta);
        dDdsigma += (bgk_stddev * bgk_stddev) / (R_cx * R_cx) *
                    (1 - std::exp(-R_cx_theta));
        diffusion_induced_drift *= dDdsigma;
    }
    TriangleMesh* triangle_mesh =
        dynamic_cast<TriangleMesh*>(const_cast<Mesh*>(&mesh));
    if (triangle_mesh != nullptr) {
        Vec3 center0 = triangle_mesh->get_center(cell_id);

        for (index_t i = 0; i < 3; i++) {
            index_t neighbour_id = triangle_mesh->get_neighbour(cell_id, i);
            if (neighbour_id == -1) {
                continue;  // No neighbour, skip
            }

            real_t R_cx_neighbour = bgk_processes[0]->rate[neighbour_id];

            Vec3 normal =
                triangle_mesh->get_surface_normal(cell_id, i, center0);
            Vec3 neighbour_center = triangle_mesh->get_center(neighbour_id);
            real_t distance = (neighbour_center - center0).norm();
            diffusion_induced_drift +=
                (R_cx_neighbour - R_cx) / distance * normal;
        }

        real_t dDdsigma = (-2 * bgk_stddev * bgk_stddev) /
                          (R_cx * R_cx * R_cx_theta) *
                          (std::exp(-R_cx_theta) - 1 + R_cx_theta);
        dDdsigma += (bgk_stddev * bgk_stddev) / (R_cx * R_cx) *
                    (1 - std::exp(-R_cx_theta));
        diffusion_induced_drift *= dDdsigma;
    }

    return (u + diffusion_induced_drift);
}

real_t calculate_diffusion_coef(
    const Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    real_t dt,
    index_t cell_id)
{
    Vec3 u = Vec3(bgk_processes[0]->bgk_mean_x[cell_id],
                  bgk_processes[0]->bgk_mean_y[cell_id],
                  bgk_processes[0]->bgk_mean_z[cell_id]);

    real_t R_cx = bgk_processes[0]->rate[cell_id];
    real_t R_cx_theta = R_cx * dt;

    real_t bgk_stddev = bgk_processes[0]->bgk_stddev[cell_id];

    real_t D = bgk_stddev * bgk_stddev / (R_cx * R_cx_theta);
    D *= (std::exp(-R_cx_theta) - 1 + R_cx_theta);
    // D *= (2. * (std::exp(-R_cx_theta) - 1) +
    //       R_cx_theta * (std::exp(-R_cx_theta) + 1));
    if (D < 0) {
        D = 0;
    }
    if (R_cx == 0) D = 0;
    return D;
}

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

    // Intermediate position to determine the drift and diffusion terms
    Vec3 drift_old = calculate_drift(mesh, bgk_processes, dt, cell_id);
    Vec3 intermediate_pos = p.geometry.position + drift_old * dt / 2.0;

    index_t intermediate_cell_id = mesh.get_cell_index(intermediate_pos);
    if (intermediate_cell_id == -1) {
        // The particle has left the mesh, just use the old position for the
        // parameters
        return false;
    }

    Vec3 drift = calculate_drift(mesh, bgk_processes, dt, intermediate_cell_id);

    real_t D0 = calculate_diffusion_coef(mesh, bgk_processes, dt, cell_id);
    real_t D = std::sqrt(2 * D0 * dt);

    Vec3 diffusion_term;
    diffusion_term.x = D * r_normal(p.rng);
    diffusion_term.y = D * r_normal(p.rng);
    diffusion_term.z = D * r_normal(p.rng);

    Vec3 new_pos = p.geometry.position + drift * dt + diffusion_term;

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
        Vec3 temp_min = p.geometry.position + drift * dt - 2 * temp_diff;
        Vec3 temp_max = p.geometry.position + drift * dt + 2 * temp_diff;

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

    // Intermediate position to determine the drift and diffusion terms
    Vec3 drift_old = calculate_drift(mesh, bgk_processes, dt, cell_id);
    Vec3 intermediate_pos = p.geometry.position + drift_old * dt / 2.0;

    index_t intermediate_cell_id = mesh.get_cell_index(intermediate_pos);
    if (intermediate_cell_id == -1) {
        // The particle has left the mesh, just use the old position for the
        // parameters
        intermediate_cell_id = cell_id;
    }

    Vec3 drift = calculate_drift(mesh, bgk_processes, dt, intermediate_cell_id);

    real_t D0 = calculate_diffusion_coef(mesh, bgk_processes, dt, cell_id);
    real_t D = std::sqrt(2 * D0 * dt);

    Vec3 diffusion_term;
    diffusion_term.x = D * r_normal(p.rng);
    diffusion_term.y = D * r_normal(p.rng);
    diffusion_term.z = D * r_normal(p.rng);

    Vec3 old_pos = p.geometry.position;
    Vec3 new_pos = old_pos + drift * dt + diffusion_term;

    // Pretend like the particle is moving in a straight line from the old the
    // the new pos
    Vec3 assumed_velocity = (1 / dt) * (new_pos - p.geometry.position);
    p.geometry.velocity = assumed_velocity;

    real_t time_left = dt;
    index_t i_step = 0;
    while (true) {
        i_step++;
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

                Vec3 normal = mesh.get_surface_normal(
                    p.geometry.cell_id, crossing.surface_id, crossing.position);

                // Make sure the normal points outwards
                if (dot(normal, p.geometry.velocity) < 0) {
                    // The normal points inwards, we need to flip it
                    normal = -normal;
                }

                real_t drift_normal = dot(drift, normal);
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
                        Vec3 new_center = mirror_pos + drift * dt;

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

                        Vec3 old_center = old_pos + drift * dt;
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
                        // use accept-reject

                        Vec3 mirror_pos =
                            old_pos -
                            2 * dot(old_pos - crossing.position, normal) *
                                normal;
                        Vec3 old_center = old_pos + drift * dt;
                        Vec3 new_center = mirror_pos + drift * dt;

                        // Coefficient for p2
                        real_t A =
                            std::exp(drift_normal / D0 *
                                     dot(crossing.position - old_pos, normal)) /
                            std::sqrt(4 * M_PI * D0 * dt);

                        // Coefficient for p3
                        real_t alpha = dot(drift, normal);
                        real_t beta = -D0;
                        real_t A2 =
                            A * 2 * (alpha / beta + drift_normal / (2 * D0));
                        // real_t A2 = 2 * (alpha / beta + drift_normal / (2 *
                        // D0));

                        real_t dist_to_crossing =
                            -dot(crossing.position - new_center, normal);
                        real_t n1;
                        while (true) {
                            // Sample from p2
                            n1 = -truncated_normal_right_tail(
                                dist_to_crossing / D, p.rng);

                            real_t n1d = n1 * D;
                            // Determine p2(x) and p3(x) for the sampled n1
                            // Scaling factors are ignored because we are only
                            // interested in the ratio
                            real_t p2 =
                                A * std::exp(-(n1d * n1d) / (4 * D0 * dt));

                            // Use numerical integration to approximate p3
                            // Note: this could definitely be improved with even
                            // something like a trapezoidal rule
                            real_t p3 = 0;
                            real_t dxi = 0.1 * std::abs(beta / alpha);
                            for (index_t i = 0; i < 10; i++) {
                                real_t xi = i * dxi;
                                real_t nn1 = n1d - xi;
                                p3 += A2 *
                                      std::exp(-(nn1 * nn1) / (4 * D0 * dt)) *
                                      std::exp(alpha / beta * xi) * dxi;
                            }

                            real_t u = r_unif(p.rng);
                            if (p2 + p3 < 0) {
                                // This is never going to get accepted
                                throw std::runtime_error(
                                    "p2 + p3 < 0, this should not happen");
                            }
                            if (u < (p2 + p3) / p2) {
                                // Accept the sample
                                break;
                            }
                        }

                        // Get the two vectors that are orthogonal to the normal
                        // and to each other using Gram-Schmidt
                        Vec3 v2(r_normal(p.rng), r_normal(p.rng),
                                r_normal(p.rng));
                        v2 -= dot(v2, normal) * normal;
                        v2 /= v2.norm();
                        Vec3 v3 = normal.cross(v2);

                        real_t n2 = r_normal(p.rng);
                        real_t n3 = r_normal(p.rng);

                        Vec3 new_pos2 = new_center + n1 * D * normal +
                                        n2 * D * v2 + n3 * D * v3;

                        // This is an extra check that should not really be
                        // necessary, check if the new position is actually
                        // inside the domain
                        if (dot(new_pos2 - crossing.position, normal) >= 0) {
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
            real_t bgk_stddev = bgk_processes[0]->bgk_stddev[cell_id];
            // Draw a new velocity
            std::normal_distribution<real_t> dist_vx(
                bgk_processes[0]->bgk_mean_x[cell_id], bgk_stddev);
            std::normal_distribution<real_t> dist_vy(
                bgk_processes[0]->bgk_mean_y[cell_id], bgk_stddev);
            std::normal_distribution<real_t> dist_vz(
                bgk_processes[0]->bgk_mean_z[cell_id], bgk_stddev);
            p.geometry.velocity =
                Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));

            std::exponential_distribution<real_t> u_dist(1);
            real_t u = u_dist(p.rng);
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
    FluidEstimator& fluid_estimator,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats,
    real_t dt,
    bool use_new_version)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t num_crossings = 0;
    index_t num_collisions = 0;
    index_t num_kd_steps = 0;
    index_t num_kd_to_kinetic_switches = 0;

    // History based simulation (because its easy)
#pragma omp parallel reduction(+ : num_crossings, num_collisions,            \
                                   num_kd_steps, num_kd_to_kinetic_switches) \
    reduction(EstimatorReduction : estimator)                                \
    reduction(FluidEstimatorReduction : fluid_estimator)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp for
        for (int i_particle = 0;
             i_particle < sources[i_source]->num_particles(); i_particle++) {
            // Generate particle
            Particle particle =
                sources[i_source]->generate_particle(i_particle, mesh);

            particle.intensity = sources[i_source]->intensity() /
                                 sources[i_source]->num_particles();

            for (index_t i_step = 0; i_step < params.max_steps; i_step++) {
                num_kd_steps++;
                Particle particle_old = particle;

                // Simulate particle kinetically until a charge exchange
                // collision occurs

                // real_t t_max = particle.t + dt;
                // if (params.t_final > 0) {
                //     t_max = std::min(t_max, params.t_final);
                // }

                while (true) {
                    particle_old = particle;
                    EventType event = update_particle(
                        mesh, bgk_processes, ion_processes, boundary_conditions,
                        particle, params.t_final);

                    estimator.add_track_length(particle, particle_old, mesh);
                    if (std::abs(particle.weight) < 1.0e-4) {
                        break;
                    }
                    if (event == EventType::FINAL_TIME) {
                        break;
                    }
                    else if (event == EventType::CHARGE_EXCHANGE) {
                        num_collisions++;
                        break;
                    }
                    else if (event == EventType::CELL_CROSSING) {
                        num_crossings++;
                    }
                }

                if (std::abs(particle.weight) < 1.0e-4) {
                    break;
                }

                // Now simulate the particle fluidly for the rest of the
                // time step

                real_t t_kdmc_end = std::ceil(particle.t / dt) * dt;
                if (params.t_final > 0) {
                    t_kdmc_end = std::min(t_kdmc_end, params.t_final);
                }

                if (particle.t < t_kdmc_end) {
                    // Now determine the final free flight time

                    real_t final_free_flight_time;

                    if (use_new_version) {
                        // number of collisions is a poisson
                        // distribution
                        std::poisson_distribution<index_t> poisson(
                            bgk_processes[0]->rate[particle.geometry.cell_id] *
                            (t_kdmc_end - particle.t));
                        index_t n_collisions = poisson(particle.rng);
                        // now that we have the number of collisions,
                        // the time of the last event is a beta
                        // distribution C++ does not have a beta
                        // distribution, so we use the two gamma
                        // distribution trick
                        std::gamma_distribution<real_t> gamma1(n_collisions, 1);
                        std::gamma_distribution<real_t> gamma2(1, 1);
                        real_t X = gamma1(particle.rng);
                        real_t Y = gamma2(particle.rng);
                        final_free_flight_time =
                            (1 - X / (X + Y)) * (t_kdmc_end - particle.t);
                    }
                    else {
                        real_t R_cx =
                            bgk_processes[0]->rate[particle.geometry.cell_id];
                        real_t eps = r_unif(particle.rng);
                        final_free_flight_time = -std::log(eps) / R_cx;
                    }

                    if (final_free_flight_time > t_kdmc_end - particle.t) {
                        final_free_flight_time = t_kdmc_end - particle.t;
                    }
                    else {
                        if (use_new_version) {
                            particle_old = particle;
                            kdmc_fluid_step(mesh, bgk_processes, ion_processes,
                                            boundary_conditions, particle,
                                            t_kdmc_end - particle.t -
                                                final_free_flight_time);

                            fluid_estimator.add_fluid_step(particle,
                                                           particle_old, mesh);
                        }
                        else {
                            bool step_succesfull = kdmc_fluid_step_old(
                                mesh, bgk_processes, ion_processes,
                                boundary_conditions, particle,
                                t_kdmc_end - particle.t -
                                    final_free_flight_time);

                            if (!step_succesfull) {
                                num_kd_to_kinetic_switches++;
                                // Particle left the mesh, switch to
                                // kinetic simulation
                                final_free_flight_time =
                                    t_kdmc_end - particle.t;
                                while (true) {
                                    particle_old = particle;
                                    EventType event = update_particle(
                                        mesh, bgk_processes, ion_processes,
                                        boundary_conditions, particle,
                                        t_kdmc_end);

                                    estimator.add_track_length(
                                        particle, particle_old, mesh);
                                    if (std::abs(particle.weight) < 1.0e-4) {
                                        break;
                                    }
                                    if (event == EventType::FINAL_TIME) {
                                        break;
                                    }
                                    else if (event ==
                                             EventType::CHARGE_EXCHANGE) {
                                        num_collisions++;
                                    }
                                    else if (event ==
                                             EventType::CELL_CROSSING) {
                                        num_crossings++;
                                    }
                                }
                            }
                            else {
                                fluid_estimator.add_fluid_step(
                                    particle, particle_old, mesh);
                            }
                        }
                    }
                    if (std::abs(particle.weight) < 1.0e-4) {
                        break;
                    }

                    // This is commented out, because instead of actually taking
                    // the free flight time as fully free flight (which would
                    // ignore the heterogeneity of the collision rate)
                    // we can just not do anything, not update the time, and the
                    // next kinetic step will take care of the kinetic part.

                    // Final free flight time
                    while (particle.t < t_kdmc_end) {
                        particle_old = particle;
                        // Create empty collision lists to simulate no
                        // collisions
                        vector<std::unique_ptr<BGKScatteringProcess>>
                            empty_bgk_processes;
                        // vector<std::unique_ptr<IonisationProcess>>
                        //     empty_ion_processes;
                        EventType event = update_particle(
                            mesh, empty_bgk_processes, ion_processes,
                            boundary_conditions, particle, t_kdmc_end);

                        estimator.add_track_length(particle, particle_old,
                                                   mesh);
                        if (std::abs(particle.weight) < 1.0e-4) {
                            break;
                        }
                        if (event == EventType::FINAL_TIME) {
                            break;
                        }
                        else if (event == EventType::CHARGE_EXCHANGE) {
                            num_collisions++;
                        }
                        else if (event == EventType::CELL_CROSSING) {
                            num_crossings++;
                        }
                    }
                }

                if (std::abs(particle.weight) < 1.0e-4) {
                    break;
                }

                // Tally
                if (params.t_final > 0) {
                    if (particle.t >= params.t_final) {
                        estimator.add_final(particle, mesh);
                        break;
                    }
                }
            }
        }
    }
}

}  // namespace neptune