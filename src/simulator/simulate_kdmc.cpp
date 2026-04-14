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
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    real_t dt,
    index_t cell_id)
{
    Vec3 u = Vec3(bgk_processes[0]->bgk_mean_x[cell_id],
                  bgk_processes[0]->bgk_mean_y[cell_id],
                  bgk_processes[0]->bgk_mean_z[cell_id]);

    real_t R_t = bgk_processes[0]->rate[cell_id];
    if (ion_processes.size() > 0) {
        R_t += ion_processes[0]->rate[cell_id];
    }
    real_t R_t_theta = R_t * dt;

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

        if (ion_processes.size() > 0) {
            if (ij.i > 0 and ij.i < voxel_mesh->nx - 1) {
                diffusion_induced_drift.x +=
                    (ion_processes[0]->rate[voxel_mesh->tuple_to_linear(
                         Index3D(ij.i + 1, ij.j, ij.k))] -
                     ion_processes[0]->rate[voxel_mesh->tuple_to_linear(
                         Index3D(ij.i - 1, ij.j, ij.k))]) /
                    (2 * voxel_mesh->dx);
            }

            if (ij.j > 0 and ij.j < voxel_mesh->ny - 1) {
                diffusion_induced_drift.y +=
                    (ion_processes[0]->rate[voxel_mesh->tuple_to_linear(
                         Index3D(ij.i, ij.j + 1, ij.k))] -
                     ion_processes[0]->rate[voxel_mesh->tuple_to_linear(
                         Index3D(ij.i, ij.j - 1, ij.k))]) /
                    (2 * voxel_mesh->dy);
            }

            if (ij.k > 0 and ij.k < voxel_mesh->nz - 1) {
                diffusion_induced_drift.z +=
                    (ion_processes[0]->rate[voxel_mesh->tuple_to_linear(
                         Index3D(ij.i, ij.j, ij.k + 1))] -
                     ion_processes[0]->rate[voxel_mesh->tuple_to_linear(
                         Index3D(ij.i, ij.j, ij.k - 1))]) /
                    (2 * voxel_mesh->dz);
            }
        }

        real_t dDdsigma = (-2 * bgk_stddev * bgk_stddev) /
                          (R_t * R_t * R_t_theta) *
                          (std::exp(-R_t_theta) - 1 + R_t_theta);
        dDdsigma += (bgk_stddev * bgk_stddev) / (R_t * R_t) *
                    (1 - std::exp(-R_t_theta));
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

            real_t R_t_neighbour = bgk_processes[0]->rate[neighbour_id];
            if (ion_processes.size() > 0) {
                R_t_neighbour += ion_processes[0]->rate[neighbour_id];
            }

            Vec3 normal =
                triangle_mesh->get_surface_normal(cell_id, i, center0);
            Vec3 neighbour_center = triangle_mesh->get_center(neighbour_id);
            real_t distance = (neighbour_center - center0).norm();
            diffusion_induced_drift +=
                (R_t_neighbour - R_t) / distance * normal;
        }

        real_t dDdsigma = (-2 * bgk_stddev * bgk_stddev) /
                          (R_t * R_t * R_t_theta) *
                          (std::exp(-R_t_theta) - 1 + R_t_theta);
        dDdsigma += (bgk_stddev * bgk_stddev) / (R_t * R_t) *
                    (1 - std::exp(-R_t_theta));
        diffusion_induced_drift *= dDdsigma;
    }

    return (u + diffusion_induced_drift);
}

real_t calculate_diffusion_coef(
    const Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    real_t dt,
    index_t cell_id)
{
    Vec3 u = Vec3(bgk_processes[0]->bgk_mean_x[cell_id],
                  bgk_processes[0]->bgk_mean_y[cell_id],
                  bgk_processes[0]->bgk_mean_z[cell_id]);

    real_t R_t = bgk_processes[0]->rate[cell_id];
    if (ion_processes.size() > 0) {
        R_t += ion_processes[0]->rate[cell_id];
    }
    real_t R_t_theta = R_t * dt;

    real_t bgk_stddev = bgk_processes[0]->bgk_stddev[cell_id];

    real_t D = bgk_stddev * bgk_stddev / (R_t * R_t_theta);
    D *= (std::exp(-R_t_theta) - 1 + R_t_theta);
    // D *= (2. * (std::exp(-R_t_theta) - 1) +
    //       R_t_theta * (std::exp(-R_t_theta) + 1));
    if (D < 0) {
        D = 0;
    }
    if (R_t == 0) D = 0;
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
    Vec3 drift_old =
        calculate_drift(mesh, bgk_processes, ion_processes, dt, cell_id);
    Vec3 intermediate_pos = p.geometry.position + drift_old * dt / 2.0;

    index_t intermediate_cell_id = mesh.get_cell_index(intermediate_pos);
    if (intermediate_cell_id == -1) {
        // The particle has left the mesh, just use the old position for the
        // parameters
        return false;
    }

    Vec3 drift = calculate_drift(mesh, bgk_processes, ion_processes, dt,
                                 intermediate_cell_id);

    real_t D0 = calculate_diffusion_coef(mesh, bgk_processes, ion_processes, dt,
                                         cell_id);
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
    Vec3 drift_old =
        calculate_drift(mesh, bgk_processes, ion_processes, dt, cell_id);
    Vec3 intermediate_pos = p.geometry.position + drift_old * dt / 2.0;

    index_t intermediate_cell_id = mesh.get_cell_index(intermediate_pos);
    if (intermediate_cell_id == -1) {
        // The particle has left the mesh, just use the old position for the
        // parameters
        intermediate_cell_id = cell_id;
    }

    Vec3 drift = calculate_drift(mesh, bgk_processes, ion_processes, dt,
                                 intermediate_cell_id);

    real_t D0 = calculate_diffusion_coef(mesh, bgk_processes, ion_processes, dt,
                                         cell_id);
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
                // Alpha and beta are the coefficients of the robin bc:
                // alpha * p + beta * dp/dn = 0
                real_t alpha, beta;
                AbsorbingBC* absorbing_bc = dynamic_cast<AbsorbingBC*>(
                    boundary_conditions[crossing.boundary_condition_id].get());
                if (absorbing_bc != nullptr) {
                    //
                    // Calculate robin BC coefficients based on maxwellian
                    // approximation of velocity distribution at the boundary.
                    //

                    real_t boundary_bgk_stddev =
                        bgk_processes[0]->bgk_stddev[p.geometry.cell_id];
                    Vec3 boundary_drift_3d =
                        Vec3(bgk_processes[0]->bgk_mean_x[p.geometry.cell_id],
                             bgk_processes[0]->bgk_mean_y[p.geometry.cell_id],
                             bgk_processes[0]->bgk_mean_z[p.geometry.cell_id]);
                    real_t boundary_drift = dot(boundary_drift_3d, normal);

                    real_t flux =
                        boundary_bgk_stddev / std::sqrt(2 * M_PI) *
                            std::exp(-0.5 *
                                     (boundary_drift / boundary_bgk_stddev) *
                                     (boundary_drift / boundary_bgk_stddev)) +
                        boundary_drift * 0.5 *
                            (1 +
                             std::erf(boundary_drift /
                                      (boundary_bgk_stddev * std::sqrt(2))));
                    alpha = drift_normal - flux;
                    beta = -D0;
                    // alpha = 1;
                    // beta = 0;
                }
                else {
                    SpecularBC* reflecting_bc = dynamic_cast<SpecularBC*>(
                        boundary_conditions[crossing.boundary_condition_id]
                            .get());
                    if (reflecting_bc != nullptr) {
                        alpha = drift_normal;
                        beta = -D0;
                    }
                    else {
                        PartiallySpecularBC* partially_reflecting_bc =
                            dynamic_cast<PartiallySpecularBC*>(
                                boundary_conditions[crossing
                                                        .boundary_condition_id]
                                    .get());
                        if (partially_reflecting_bc != nullptr) {
                            // For a partially reflecting BC, we do something in between
                            // specular and absorbing
                            real_t boundary_bgk_stddev =
                                bgk_processes[0]->bgk_stddev[p.geometry.cell_id];
                            Vec3 boundary_drift_3d =
                                Vec3(bgk_processes[0]->bgk_mean_x[p.geometry.cell_id],
                                    bgk_processes[0]->bgk_mean_y[p.geometry.cell_id],
                                    bgk_processes[0]->bgk_mean_z[p.geometry.cell_id]);
                            real_t boundary_drift = dot(boundary_drift_3d, normal);

                            real_t flux =
                                boundary_bgk_stddev / std::sqrt(2 * M_PI) *
                                    std::exp(-0.5 *
                                            (boundary_drift / boundary_bgk_stddev) *
                                            (boundary_drift / boundary_bgk_stddev)) +
                                boundary_drift * 0.5 *
                                    (1 +
                                    std::erf(boundary_drift /
                                            (boundary_bgk_stddev * std::sqrt(2))));

                            alpha = drift_normal - flux * partially_reflecting_bc->absorption_coefficient;
                            beta = -D0;
                        }
                        else {
                            throw std::runtime_error(
                                "The chosen boundary condition is not "
                                "supported "
                                "for kdmc");
                        }
                    }
                }

                Vec3 mirror_pos =
                    old_pos -
                    2 * dot(old_pos - crossing.position, normal) * normal;
                Vec3 old_center = old_pos + drift * dt;
                Vec3 new_center = mirror_pos + drift * dt;

                // Coefficient for p2
                real_t A = std::exp(drift_normal / D0 *
                                    dot(crossing.position - old_pos, normal));

                // Coefficient for p3
                real_t A2;
                if (beta == 0) {
                    A2 = 0;
                    A = -A;
                }
                else {
                    A2 = -2 * A * (alpha / beta + drift_normal / (2 * D0));
                }

                real_t dist_to_crossing =
                    -dot(crossing.position - new_center, normal);
                // Sample from p2
                real_t n1 =
                    -truncated_normal_right_tail(dist_to_crossing / D, p.rng);

                real_t n1d = n1 * D;
                real_t s = new_center.x + n1d * normal.x;
                // Determine p2(x) and p3(x) for the sampled n1
                real_t p2 = A / std::sqrt(4 * M_PI * D0 * dt) *
                            std::exp(-(n1d * n1d) / (4 * D0 * dt));

                real_t p3 = 0;
                if (beta != 0 and std::abs(A2) > 1.0e-8) {
                    real_t sqrt4Dt = std::sqrt(4.0 * D0 * dt);
                    real_t z =
                        (-n1d + 2.0 * (alpha / beta) * D0 * dt) / sqrt4Dt;
                    real_t B = (alpha / beta) * (-n1d) +
                               (alpha / beta) * (alpha / beta) * D0 * dt;
                    if (B > 700) {
                        // To prevent overflow in the exp, we can use the
                        // asymptotic expansion of the erfc for large z
                        p3 = A2 * 0.5 * std::exp(B - z * z) /
                             (z * std::sqrt(M_PI));
                    }
                    else {
                        p3 = A2 * 0.5 * std::exp(B) * std::erfc(z);
                    }
                    // Check if p3 is NaN, if so, throw an error
                    if (std::isnan(p3)) {
                        throw std::runtime_error(
                            "p3 is NaN, this should not happen, check the "
                            "calculations for p3");
                    }
                }

                // Multiply particle weight by int_-inf^L p2(x)dx / int_L^inf
                // p1(x)dx to account for the area's under p1 and p2
                real_t dist_to_crossing_p1 =
                    -dot(crossing.position - old_center, normal);
                real_t area_p1 =
                    0.5 * std::erfc(-dist_to_crossing_p1 / (D * std::sqrt(2)));
                real_t area_p2 = A * (1 - 0.5 * std::erfc(-dist_to_crossing /
                                                          (D * std::sqrt(2))));
                p.weight *= (area_p2) / (area_p1);
                // Multiply particle weight by (p2 + p3) / p2 to account for the
                // fact that we are only sampling from p2 for the reflected
                // particles
                p.weight *= (p2 + p3) / p2;

                // Get the two vectors that are orthogonal to the normal
                // and to each other using Gram-Schmidt
                Vec3 v2(r_normal(p.rng), r_normal(p.rng), r_normal(p.rng));
                v2 -= dot(v2, normal) * normal;
                v2 /= v2.norm();
                Vec3 v3 = normal.cross(v2);

                real_t n2 = r_normal(p.rng);
                real_t n3 = r_normal(p.rng);

                Vec3 new_pos2 =
                    new_center + n1 * D * normal + n2 * D * v2 + n3 * D * v3;

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
    stats.num_crossings = num_crossings;
    stats.num_collisions = num_collisions;
    stats.num_kd_steps = num_kd_steps;
    stats.num_kd_to_kinetic_switches = num_kd_to_kinetic_switches;
    stats.total_particles = 0.0;
    for (const auto& source : sources) {
        stats.total_particles += source->num_particles();
    }
}

}  // namespace neptune