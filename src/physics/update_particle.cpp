#include "neptune/physics/update_particle.hpp"

#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

EventType update_particle(
    const Mesh& mesh,
    const PlasmaBackground& background,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    Particle& p,
    real_t t_final)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    // Calculate total collision rate (charge exchange + neutral-neutral)

    real_t R_cx = constants::calculate_K_cx(background.plasma_temp[cell_id]) *
                  background.plasma_dens[cell_id];
    real_t R_nn = constants::calculate_K_cx(background.neutral_temp[cell_id]) *
                  background.neutral_dens[cell_id];

    real_t collision_rate = R_cx + R_nn;

    // random number to determine time till collision
    real_t eps = r_unif(p.rng);
    // time till collision
    real_t tc = -std::log(eps) / collision_rate;

    // Calculate the next cell crossing
    CellCrossing crossing = mesh.get_next_crossing(p.geometry);

    // Check which one is earlier, charge exchange, cell crossing or final time
    if (t_final > 0) {
        real_t time_left = t_final - p.t;
        if (time_left < tc && time_left < crossing.t) {
            // We will not hit the boundary before the time is up
            // This is the final iteration
            p.geometry.position += time_left * p.geometry.velocity;
            p.geometry.boundary_id = -1;
            p.t = t_final;
            p.weight =
                p.weight * std::exp(-time_left *
                                    constants::calculate_K_i(
                                        background.plasma_temp[cell_id]) *
                                    background.plasma_dens[cell_id]);
            return EventType::FINAL_TIME;
        }
    }
    real_t dt;
    EventType event_type;
    if (tc < crossing.t) {
        event_type = EventType::CHARGE_EXCHANGE;
        // Move the particle until the collision
        p.geometry.position += tc * p.geometry.velocity;

        // Determine which collision happened, charge exchange or
        // neutral-neutral
        real_t eps = r_unif(p.rng);

        if (eps < R_cx / collision_rate) {
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
        }
        else {
            // Generate new velocity from neutral background (BGK
            // operator)
            real_t stddev = constants::calculate_BGK_stddev(
                background.neutral_temp[cell_id]);
            std::normal_distribution<real_t> dist_vx(
                background.neutral_vel_x[cell_id], stddev);
            std::normal_distribution<real_t> dist_vy(
                background.neutral_vel_y[cell_id], stddev);
            std::normal_distribution<real_t> dist_vz(
                background.neutral_vel_z[cell_id], stddev);
            p.geometry.velocity =
                Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));
        }

        // Reset boundary id (we are no longer at a boundary)
        p.geometry.boundary_id = -1;

        dt = tc;
    }
    else {
        event_type = EventType::CELL_CROSSING;
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

            Vec3 normal = mesh.get_surface_normal(cell_id, crossing.surface_id);

            boundary_conditions[crossing.boundary_condition_id]->apply(p,
                                                                       normal);
        }
    }

    // Ionisation
    p.weight =
        p.weight *
        std::exp(-dt *
                 constants::calculate_K_i(background.plasma_temp[cell_id]) *
                 background.plasma_dens[cell_id]);

    p.t += dt;

    return event_type;
}

}  // namespace neptune