#include <omp.h>

#include "eigrene/geometry/voxelmesh.hpp"
#include "eigrene/physics/linesource.hpp"
#include "eigrene/physics/physicalconstants.hpp"
#include "eigrene/physics/update_particle.hpp"
#include "eigrene/simulator/simulator.hpp"

namespace eigrene {

/**
 * The kinetic part of a kdmc update
 *
 * return true if the last event was a collision
 * return false if the last event was a crossing/boundary event
 */
bool update_particle_kdmc_kinetic(
    const Mesh& mesh,
    const PlasmaBackground& background,
    const vector<BoundaryCondition>& boundary_conditions,
    Particle& p,
    index_t& num_crossings,
    index_t& num_collisions,
    real_t t_final)
{
    bool ended_in_collision = false;
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    // random number to determine time till collision
    real_t eps = r_unif(p.rng);
    // time till collision
    real_t tc = -std::log(eps) /
                (constants::calculate_K_cx(background.plasma_temp[cell_id]) *
                 background.plasma_dens[cell_id]);

    // Calculate the next cell crossing
    CellCrossing crossing = mesh.get_next_crossing(p.geometry);

    // Check which one is earlier, charge exchange, cell crossing or final time
    real_t time_left = t_final - p.t;
    if (time_left < tc && time_left < crossing.t) {
        // We will not hit the boundary before the time is up
        p.geometry.position += time_left * p.geometry.velocity;
        p.geometry.boundary_id = -1;
        p.t = t_final;
        p.weight =
            p.weight *
            std::exp(-time_left *
                     constants::calculate_K_i(background.plasma_temp[cell_id]) *
                     background.plasma_dens[cell_id]);
        return;
    }
    real_t dt;
    if (tc < crossing.t) {
        num_collisions++;
        ended_in_collision = true;
        // Move the particle until the collision
        p.geometry.position += tc * p.geometry.velocity;
        // Generate new velocity from plasma background (BGK
        // operator)
        real_t stddev =
            constants::calculate_BGK_stddev(background.plasma_temp[cell_id]);
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

        dt = tc;
    }
    else {
        num_crossings++;
        dt = crossing.t;
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
            }
            else {
                throw std::runtime_error(
                    "Using a boundary condition that is not (yet) implemented");
            }
        }
    }

    p.t += dt;

    return ended_in_collision;
}

real_t sample_inverse_gaussian(real_t mu, real_t lambda, rng_t& rng)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);
    real_t n = r_normal(rng);
    real_t y = n * n;
    real_t x =
        mu + (mu * mu * y) / (2 * lambda) -
        (mu / (2 * lambda)) * std::sqrt(4 * mu * lambda * y + mu * mu * y * y);

    real_t u = r_unif(rng);

    if (u <= mu / (mu + x)) {
        return x;
    }
    else {
        return mu * mu / x;
    }
}

std::pair<real_t, index_t> sample_first_exit_time(
    const Mesh& mesh, index_t cell_id, Vec3 pos, Vec3 u, real_t D, rng_t& rng)
{
    VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);

    real_t dx = vmesh.dx;
    real_t dy = vmesh.dy;
    real_t dz = vmesh.dz;

    Index3D ij = vmesh.linear_to_tuple(cell_id);

    index_t boundary = -1;
    real_t t = std::numeric_limits<real_t>::infinity();

    // check x 1
    real_t d = pos.x - ij.i * dx;
    real_t v = -u.x;
    real_t t_temp;
    if (std::abs(v) > 0) {
        t_temp = sample_inverse_gaussian(d / v, (d * d) / (2 * D), rng);
    }
    else {
        t_temp = sample_inverse_gaussian(d / D, (d * d) / (2 * D), rng);
    }
    if (t_temp < t) {
        t = t_temp;
        boundary = 0;
    }
    // check y 1
    d = pos.y - ij.j * dy;
    v = -u.y;
    if (std::abs(v) > 0) {
        t_temp = sample_inverse_gaussian(d / v, (d * d) / (2 * D), rng);
    }
    else {
        t_temp = sample_inverse_gaussian(d / D, (d * d) / (2 * D), rng);
    }
    if (t_temp < t) {
        t = t_temp;
        boundary = 1;
    }
    // check z 1
    d = pos.z - ij.k * dz;
    v = -u.z;
    if (std::abs(v) > 0) {
        t_temp = sample_inverse_gaussian(d / v, (d * d) / (2 * D), rng);
    }
    else {
        t_temp = sample_inverse_gaussian(d / D, (d * d) / (2 * D), rng);
    }
    if (t_temp < t) {
        t = t_temp;
        boundary = 2;
    }
    // check x 2
    d = (ij.i + 1) * dx - pos.x;
    v = u.x;
    if (std::abs(v) > 0) {
        t_temp = sample_inverse_gaussian(d / v, (d * d) / (2 * D), rng);
    }
    else {
        t_temp = sample_inverse_gaussian(d / D, (d * d) / (2 * D), rng);
    }
    if (t_temp < t) {
        t = t_temp;
        boundary = 3;
    }
    // check y 2
    d = (ij.j + 1) * dy - pos.y;
    v = u.y;
    if (std::abs(v) > 0) {
        t_temp = sample_inverse_gaussian(d / v, (d * d) / (2 * D), rng);
    }
    else {
        t_temp = sample_inverse_gaussian(d / D, (d * d) / (2 * D), rng);
    }
    if (t_temp < t) {
        t = t_temp;
        boundary = 4;
    }
    // check z 2
    d = (ij.k + 1) * dz - pos.z;
    v = u.z;
    if (std::abs(v) > 0) {
        t_temp = sample_inverse_gaussian(d / v, (d * d) / (2 * D), rng);
    }
    else {
        t_temp = sample_inverse_gaussian(d / D, (d * d) / (2 * D), rng);
    }
    if (t_temp < t) {
        t = t_temp;
        boundary = 5;
    }

    return std::make_pair(t, boundary);
}

void update_particle_kdmc_diffusive(
    const Mesh& mesh,
    const PlasmaBackground& background,
    const vector<BoundaryCondition>& boundary_conditions,
    Particle& p,
    index_t& num_crossings,
    index_t& num_collisions,
    real_t t_final,
    real_t dt,
    real_t t_kinetic)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    real_t diffusive_flight_time = dt - t_kinetic;
    diffusive_flight_time = std::min(diffusive_flight_time, t_final - p.t);

    if (diffusive_flight_time <= 0) {
        return;
    }

    //
    // Calculate the drift and diffusion terms.
    // This is equation 3.1 -- 3.6 in the KDMC paper.
    // We don't do scaling with epsilon here, it should be safe
    // to just remove epsilon from the equations.
    //

    // Charge exchange rate coefficient at the intermediate cell
    real_t K_cx = constants::calculate_K_cx(background.plasma_temp[cell_id]);
    real_t R_cx = K_cx * background.plasma_dens[cell_id];
    real_t R_cx_theta = R_cx * diffusive_flight_time;

    Vec3 u =
        Vec3(background.plasma_vel_x[cell_id], background.plasma_vel_y[cell_id],
             background.plasma_vel_z[cell_id]);

    real_t D = background.plasma_temp[cell_id] * constants::eV_to_Joules /
               constants::m_p / (R_cx * R_cx_theta);
    D *= (std::exp(-R_cx_theta) - 1 + R_cx_theta);

    if (D < 0) {
        D = 0;
    }

    // Calculate times to hit the boundaries
    // This only works for voxelmeshes for now although the math should also
    // work for triangles
    auto [t_exit, boundary] = sample_first_exit_time(
        mesh, cell_id, p.geometry.position, u, D, p.rng);

    // Check if we hit a boundary before the time is up
    if (diffusive_flight_time < t_exit) {
        // draw a new position, but make sure it is within the current cell
        while(true){
            Vec3 new_pos = p.geometry.position + diffusive_flight_time * u + std::sqrt(2 * D * diffusive_flight_time) * Vec3(r_normal(p.rng), r_normal(p.rng), r_normal(p.rng));
            if(mesh.is_in_cell(new_pos, cell_id)){
                p.geometry.position = new_pos;
                break;
            }
        }
        p.t += diffusive_flight_time;
    } else {
        // Sample the position at the boundary
        // one coordinate is already fixed by the boundary
        // we assume the other are independent

        VoxelMesh vmesh = dynamic_cast<const VoxelMesh&>(mesh);
        Index3D ij = vmesh.linear_to_tuple(cell_id);
        Vec3 new_pos;
        if(boundary == 0 or boundary == 3){
            new_pos.x = (boundary == 0) ? ij.i * vmesh.dx : (ij.i + 1) * vmesh.dx;
            while(true){
                new_pos.y = p.geometry.position.y + t_exit * u.y + std::sqrt(2 * D * t_exit) * r_normal(p.rng);
                if(new_pos.y >= ij.j * vmesh.dy and new_pos.y <= (ij.j + 1) * vmesh.dy){
                    break;
                }
            }
            while(true){
                new_pos.z = p.geometry.position.z + t_exit * u.z + std::sqrt(2 * D * t_exit) * r_normal(p.rng);
                if(new_pos.z >= ij.k * vmesh.dz and new_pos.z <= (ij.k + 1) * vmesh.dz){
                    break;
                }
            }
        }
        else if(boundary == 1 or boundary == 4){
            new_pos.y = (boundary == 1) ? ij.j * vmesh.dy : (ij.j + 1) * vmesh.dy;
            while(true){
                new_pos.x = p.geometry.position.x + t_exit * u.x + std::sqrt(2 * D * t_exit) * r_normal(p.rng);
                if(new_pos.x >= ij.i * vmesh.dx and new_pos.x <= (ij.i + 1) * vmesh.dx){
                    break;
                }
            }
            while(true){
                new_pos.z = p.geometry.position.z + t_exit * u.z + std::sqrt(2 * D * t_exit) * r_normal(p.rng);
                if(new_pos.z >= ij.k * vmesh.dz and new_pos.z <= (ij.k + 1) * vmesh.dz){
                    break;
                }
            }
        }
        else if(boundary == 2 or boundary == 5){
            new_pos.z = (boundary == 2) ? ij.k * vmesh.dz : (ij.k + 1) * vmesh.dz;
            while(true){
                new_pos.x = p.geometry.position.x + t_exit * u.x + std::sqrt(2 * D * t_exit) * r_normal(p.rng);
                if(new_pos.x >= ij.i * vmesh.dx and new_pos.x <= (ij.i + 1) * vmesh.dx){
                    break;
                }
            }
            while(true){
                new_pos.y = p.geometry.position.y + t_exit * u.y + std::sqrt(2 * D * t_exit) * r_normal(p.rng);
                if(new_pos.y >= ij.j * vmesh.dy and new_pos.y <= (ij.j + 1) * vmesh.dy){
                    break;
                }
            }
        }
        p.geometry.position = new_pos;
        p.t += t_exit;
        p.geometry.boundary_id = boundary;

        // now handle boundary conditions

    }

}

void simulate_kdmc2(Mesh& mesh,
                    PlasmaBackground& background,
                    vector<std::unique_ptr<Source>>& sources,
                    vector<real_t>& density,
                    vector<BoundaryCondition>& boundary_conditions,
                    SimulationParams params,
                    ExecutionStatistics& stats,
                    real_t t_final,
                    real_t dt)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    // initialize the density array to zero
    for (index_t i = 0; i < mesh.n_cells(); i++) {
        density[i] = 0;
    }

    // History based simulation (because its easy)
    for (int i_source = 0; i_source < sources.size(); i_source++) {
#pragma omp parallel
        {
            index_t num_crossings_local = 0;
            index_t num_collisions_local = 0;

            // This should really be changed to a reduction clause instead of
            // separate local variables
            vector<real_t> density_local(mesh.n_cells(), 0);
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
                    // Simulate the particle up to the next event
                    update_particle_kdmc2(mesh, background, boundary_conditions,
                                          particle, num_crossings_local,
                                          num_collisions_local, t_final, dt);

                    // Tally
                    if (particle.t >= t_final) {
                        real_t v = mesh.get_volume(particle.geometry.cell_id);
                        density_local[particle.geometry.cell_id] +=
                            intensity_multiplier / v;
                        break;
                    }

                    if (particle.weight < 1.0e-4) {
                        break;
                    }
                }

#ifdef _OPENMP
                auto thread_id = omp_get_thread_num();
                if (thread_id == 0) {
                    if (i_particle % 1000 == 0) {
                        std::cout.flush();
                        std::cout
                            << "Progress: "
                            << i_particle /
                                   ((float)sources[i_source]->num_particles()) *
                                   100
                            << "% \r";
                    }
                }
#endif
            }

            for (index_t i = 0; i < mesh.n_cells(); i++) {
#pragma omp atomic
                density[i] += density_local[i];
            }
#pragma omp atomic
            stats.num_crossings += num_crossings_local;
#pragma omp atomic
            stats.num_collisions += num_collisions_local;
        }
    }
}

}  // namespace eigrene