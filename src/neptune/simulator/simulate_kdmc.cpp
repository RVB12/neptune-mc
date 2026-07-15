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

/**
 * Calculate the drift term used in the KDMC method.
 * This is not merely the mean post-collision velocity, but also includes drift
 * induced by slopes in the collision rate.
 * Note, this slope currently only works for voxel meshes
 */
Vec3 calculate_drift_kdmc(
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
    real_t dt,
    index_t cell_id)
{
    // Make a pretend particle with the background species and the cell index
    Particle particle(0, 0, cx_collisions[0]->particle_species());
    particle.geometry.cell_id = cell_id;

    Vec3 u =
        cx_collisions[0]->post_collision_velocity_mean(particle, background);

    real_t R_t = 0;
    for (index_t i = 0; i < cx_collisions.size(); i++) {
        R_t += cx_collisions[i]->collision_rate(particle, background);
    }
    for (index_t i = 0; i < ion_collisions.size(); i++) {
        R_t += ion_collisions[i]->collision_rate(particle, background);
    }

    real_t R_t_theta = R_t * dt;

    real_t bgk_stddev =
        cx_collisions[0]->post_collision_velocity_var(particle, background);
    bgk_stddev = std::sqrt(bgk_stddev);

    // Calculate derivative of the diffusion term and add it to the drift term
    Vec3 diffusion_induced_drift(0, 0, 0);

    Vec3 center0 = mesh.get_center(cell_id);

    for (index_t i = 0; i < mesh.get_n_surfaces(cell_id); i++) {
        index_t neighbour_id = mesh.get_neighbour(cell_id, i);
        if (neighbour_id < 0) {
            continue;  // No neighbour, skip
        }

        particle.geometry.cell_id = neighbour_id;
        real_t R_t_neighbour = 0;
        for (index_t i = 0; i < cx_collisions.size(); i++) {
            R_t_neighbour +=
                cx_collisions[i]->collision_rate(particle, background);
        }
        for (index_t i = 0; i < ion_collisions.size(); i++) {
            R_t_neighbour +=
                ion_collisions[i]->collision_rate(particle, background);
        }

        Vec3 normal = mesh.get_surface_normal(cell_id, i, center0);
        Vec3 neighbour_center = mesh.get_center(neighbour_id);
        real_t distance = (neighbour_center - center0).norm();
        diffusion_induced_drift += (R_t_neighbour - R_t) / distance * normal;
    }

    real_t dDdRt = (-2 * bgk_stddev * bgk_stddev) / (R_t * R_t * R_t_theta) *
                   (std::exp(-R_t_theta) - 1 + R_t_theta);
    dDdRt +=
        (bgk_stddev * bgk_stddev) / (R_t * R_t) * (1 - std::exp(-R_t_theta));
    diffusion_induced_drift *= dDdRt;

    return (u + diffusion_induced_drift);
}

real_t calculate_diffusion_coef(
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
    real_t dt,
    index_t cell_id)
{
    // Make a pretend particle with the background species and the cell index
    Particle particle(0, 0, cx_collisions[0]->particle_species());
    particle.geometry.cell_id = cell_id;

    Vec3 u =
        cx_collisions[0]->post_collision_velocity_mean(particle, background);

    real_t R_t = 0;
    for (index_t i = 0; i < cx_collisions.size(); i++) {
        R_t += cx_collisions[i]->collision_rate(particle, background);
    }
    for (index_t i = 0; i < ion_collisions.size(); i++) {
        R_t += ion_collisions[i]->collision_rate(particle, background);
    }
    real_t R_t_theta = R_t * dt;

    real_t bgk_stddev =
        cx_collisions[0]->post_collision_velocity_var(particle, background);
    bgk_stddev = std::sqrt(bgk_stddev);

    // Note: this formula is differs from Bert's paper, that is because
    // we don't use KDMC, but rather KDKMC and we don't need the correction
    // for the final velocity.
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
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
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
    Vec3 drift_old = calculate_drift_kdmc(mesh, background, cx_collisions,
                                          ion_collisions, dt, cell_id);
    Vec3 intermediate_pos = p.geometry.position + drift_old * dt / 2.0;

    index_t intermediate_cell_id = mesh.get_cell_index(intermediate_pos);
    if (intermediate_cell_id == -1) {
        // The particle has left the mesh, just use the old position for the
        // parameters
        intermediate_cell_id = cell_id;
        intermediate_pos = p.geometry.position;
        // return false;
    }

    Vec3 drift = calculate_drift_kdmc(mesh, background, cx_collisions,
                                      ion_collisions, dt, intermediate_cell_id);

    real_t D0 = calculate_diffusion_coef(mesh, background, cx_collisions,
                                         ion_collisions, dt, cell_id);
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
    if (ion_collisions.size() > 0) {
        p.weight *=
            std::exp(-dt * ion_collisions[0]->collision_rate(p, background));
    }
    p.geometry.position = new_pos;
    p.geometry.cell_id = new_cell_id;
    p.geometry.boundary_id = -1;
    p.t += dt;
    return true;
}

/**
 * @brief Simulate a particle using a diffusive fluid approximation
 *
 * We assume there is only one bgk process
 */
void kdmc_fluid_step(
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
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
    Vec3 drift_old = calculate_drift_kdmc(mesh, background, cx_collisions,
                                          ion_collisions, dt, cell_id);
    Vec3 intermediate_pos = p.geometry.position + drift_old * dt / 2.0;

    index_t intermediate_cell_id = mesh.get_cell_index(intermediate_pos);
    if (intermediate_cell_id == -1) {
        // The particle has left the mesh, just use the old position for the
        // parameters
        intermediate_cell_id = cell_id;
    }

    Vec3 drift = calculate_drift_kdmc(mesh, background, cx_collisions,
                                      ion_collisions, dt, intermediate_cell_id);

    real_t D0 = calculate_diffusion_coef(mesh, background, cx_collisions,
                                         ion_collisions, dt, cell_id);
    real_t D = std::sqrt(2 * D0 * dt);

    Vec3 diffusion_term;
    diffusion_term.x = D * r_normal(p.rng);
    diffusion_term.y = D * r_normal(p.rng);
    diffusion_term.z = D * r_normal(p.rng);

    Vec3 old_pos = p.geometry.position;
    Vec3 new_pos = old_pos + drift * dt + diffusion_term;

    // Pretend like the particle is moving in a straight line from the old to
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
                        cx_collisions[0]->post_collision_velocity_var(
                            p, background);
                    boundary_bgk_stddev = std::sqrt(boundary_bgk_stddev);

                    Vec3 boundary_drift_3d =
                        cx_collisions[0]->post_collision_velocity_mean(
                            p, background);

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
                            // For a partially reflecting BC, we do something in
                            // between specular and absorbing
                            real_t boundary_bgk_stddev =
                                cx_collisions[0]->post_collision_velocity_var(
                                    p, background);
                            boundary_bgk_stddev =
                                std::sqrt(boundary_bgk_stddev);

                            Vec3 boundary_drift_3d =
                                cx_collisions[0]->post_collision_velocity_mean(
                                    p, background);
                            real_t boundary_drift =
                                dot(boundary_drift_3d, normal);

                            real_t flux =
                                boundary_bgk_stddev / std::sqrt(2 * M_PI) *
                                    std::exp(
                                        -0.5 *
                                        (boundary_drift / boundary_bgk_stddev) *
                                        (boundary_drift /
                                         boundary_bgk_stddev)) +
                                boundary_drift * 0.5 *
                                    (1 + std::erf(boundary_drift /
                                                  (boundary_bgk_stddev *
                                                   std::sqrt(2))));

                            alpha = drift_normal -
                                    flux * partially_reflecting_bc
                                               ->absorption_coefficient;
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

                // TODO: we resample the two orthogonal components
                // but (for specular reflection) we should actually just keep
                // the same orthogonal components as before
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

            p.geometry.velocity =
                cx_collisions[0]->sample_post_collision_velocity(p, background);
            std::exponential_distribution<real_t> u_dist(1);
            real_t u = u_dist(p.rng);
            break;
        }
    }

    for (index_t i = 0; i < ion_collisions.size(); i++) {
        p.weight *=
            std::exp(-dt * ion_collisions[i]->collision_rate(p, background));
    }
}

void simulate_kdmc(
    Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
    vector<std::unique_ptr<Source>>& sources,
    EstimatorGroup& estimator_group,
    vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    SimulationParams params,
    ExecutionStatistics& stats)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // Calculate total number of particles
    index_t total_particles = 0;
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        total_particles += sources[i_source]->num_particles();
    }

    index_t num_crossings = 0;
    index_t num_collisions = 0;
    index_t num_kd_steps = 0;
    index_t num_kd_to_kinetic_switches = 0;

    stats.total_particles = total_particles * params.n_batches;
    vector<size_t> cumulative_particles_per_source(sources.size() + 1, 0);
    for (int i_source = 0; i_source < sources.size(); i_source++) {
        cumulative_particles_per_source[i_source + 1] =
            cumulative_particles_per_source[i_source] +
            sources[i_source]->num_particles();
    }

    real_t dt = params.kdmc_dt;

    // History based simulation (because its easy)
#pragma omp parallel reduction(+ : num_crossings, num_collisions) \
    reduction(EstimatorGroupReduction : estimator_group)
    for (index_t i_batch = 0; i_batch < params.n_batches; i_batch++) {
#pragma omp for schedule(guided, 100)
        for (index_t global_particle_index = 0;
             global_particle_index < total_particles; global_particle_index++) {
            // Find source and local particles index
            index_t i_source = 0;
            // binary search to find the source index
            index_t left = 0;
            index_t right = sources.size();
            while (left < right) {
                index_t mid = left + (right - left) / 2;
                if (cumulative_particles_per_source[mid] <=
                    global_particle_index) {
                    left = mid + 1;
                }
                else {
                    right = mid;
                }
            }
            i_source = left - 1;

            Particle particle = sources[i_source]->generate_particle(
                global_particle_index + i_batch * total_particles, params.seed,
                mesh, background);

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
                        mesh, background, cx_collisions, ion_collisions,
                        boundary_conditions, particle, params.t_final);

                    estimator_group.tally_track_length(particle, particle_old,
                                                       mesh, background,
                                                       ion_collisions);
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

                    real_t R_cx = 0;
                    for (index_t i = 0; i < cx_collisions.size(); i++) {
                        R_cx += cx_collisions[i]->collision_rate(particle,
                                                                 background);
                    }
                    if (params.algorithm == AlgorithmType::KDMC) {
                        // number of collisions is a poisson
                        // distribution

                        std::poisson_distribution<index_t> poisson(
                            R_cx * (t_kdmc_end - particle.t));
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
                        // Assume the final time is also exponentially
                        // distributed with the same rate as the other collision
                        // times
                        real_t eps = r_unif(particle.rng);
                        final_free_flight_time = -std::log(eps) / R_cx;
                    }

                    if (final_free_flight_time > t_kdmc_end - particle.t) {
                        final_free_flight_time = t_kdmc_end - particle.t;
                    }
                    else {
                        if (params.algorithm == AlgorithmType::KDMC) {
                            particle_old = particle;
                            kdmc_fluid_step(mesh, background, cx_collisions,
                                            ion_collisions, boundary_conditions,
                                            particle,
                                            t_kdmc_end - particle.t -
                                                final_free_flight_time);

                            estimator_group.tally_fluid(particle, particle_old,
                                                        mesh);
                        }
                        else {
                            particle_old = particle;
                            bool step_succesfull = kdmc_fluid_step_old(
                                mesh, background, cx_collisions, ion_collisions,
                                boundary_conditions, particle,
                                t_kdmc_end - particle.t -
                                    final_free_flight_time);

                            if (step_succesfull) {
                                estimator_group.tally_fluid(particle,
                                                            particle_old, mesh);
                            }
                            else {
                                num_kd_to_kinetic_switches++;
                                // Particle left the mesh, switch to
                                // kinetic simulation
                                final_free_flight_time = 0.0;
                                while (true) {
                                    particle_old = particle;
                                    EventType event = update_particle(
                                        mesh, background, cx_collisions,
                                        ion_collisions, boundary_conditions,
                                        particle, t_kdmc_end);

                                    estimator_group.tally_track_length(
                                        particle, particle_old, mesh,
                                        background, ion_collisions);
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
                        }
                    }
                    if (std::abs(particle.weight) < 1.0e-4) {
                        break;
                    }

                    // This can be commented out, because instead of actually
                    // taking the free flight time as fully free flight (which
                    // would ignore the heterogeneity of the collision rate)
                    // we can just not do anything, not update the time, and the
                    // next kinetic step will take care of the kinetic part.

                    // Final free flight time
                    // vector<std::unique_ptr<ChargeExchangeCollision>>
                    // cx_collisions_empty; while (particle.t < t_kdmc_end) {

                    //     particle_old = particle;
                    //     EventType event = update_particle(
                    //         mesh, background, cx_collisions_empty,
                    //         ion_collisions, boundary_conditions, particle,
                    //         t_kdmc_end);

                    //     estimator_group.tally_track_length(
                    //         particle, particle_old, mesh, background,
                    //         ion_collisions);
                    //     if (std::abs(particle.weight) < 1.0e-4) {
                    //         break;
                    //     }
                    //     if (event == EventType::FINAL_TIME) {
                    //         break;
                    //     }
                    //     else if (event == EventType::CHARGE_EXCHANGE) {
                    //         num_collisions++;
                    //     }
                    //     else if (event == EventType::CELL_CROSSING) {
                    //         num_crossings++;
                    //     }
                    // }
                }

                if (std::abs(particle.weight) < 1.0e-4) {
                    break;
                }

                // Tally
                if (params.t_final > 0) {
                    if (particle.t >= params.t_final) {
                        estimator_group.tally_final(particle, mesh);
                        particle.weight = 0;
                        break;
                    }
                }
            }
        }
        estimator_group.finalize_batch();
    }
    estimator_group.finalize();

    stats.num_crossings = num_crossings;
    stats.num_collisions = num_collisions;
    stats.num_kd_steps = num_kd_steps;
    stats.num_kd_to_kinetic_switches = num_kd_to_kinetic_switches;
}

}  // namespace neptune