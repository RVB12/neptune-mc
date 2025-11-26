#include "neptune/physics/update_particle.hpp"

#include "neptune/physics/physicalconstants.hpp"

namespace neptune {

EventType update_particle(
    const Mesh& mesh,
    const vector<std::unique_ptr<BGKScatteringProcess>>& bgk_processes,
    const vector<std::unique_ptr<IonisationProcess>>& ion_processes,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    Particle& p,
    real_t t_final)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    // Calculate total collision rate (charge exchange + neutral-neutral)

    real_t collision_rate = 0;
    for (int i = 0; i < bgk_processes.size(); i++) {
        collision_rate += bgk_processes[i]->rate[cell_id];
    }

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
            for (int i = 0; i < ion_processes.size(); i++) {
                p.weight *=
                    std::exp(-time_left * ion_processes[i]->rate[cell_id]);
            }
            return EventType::FINAL_TIME;
        }
    }
    real_t dt;
    EventType event_type;
    if (tc < crossing.t) {
        event_type = EventType::CHARGE_EXCHANGE;
        // Move the particle until the collision
        p.geometry.position += tc * p.geometry.velocity;

        // A collision has occurred, but we don't know which one
        // draw a new random number to determine that
        real_t eps = r_unif(p.rng);

        real_t cumulative_rate = 0;
        for (int i = 0; i < bgk_processes.size(); i++) {
            cumulative_rate += bgk_processes[i]->rate[cell_id];
            if (eps <= cumulative_rate / collision_rate) {
                // We have a collision with the background
                // Generate new velocity from plasma background (BGK
                // operator)
                real_t stddev = bgk_processes[i]->bgk_stddev[cell_id];
                std::normal_distribution<real_t> dist_vx(
                    bgk_processes[i]->bgk_mean_x[cell_id], stddev);
                std::normal_distribution<real_t> dist_vy(
                    bgk_processes[i]->bgk_mean_y[cell_id], stddev);
                std::normal_distribution<real_t> dist_vz(
                    bgk_processes[i]->bgk_mean_z[cell_id], stddev);
                Vec3 new_velocity_notransform =
                    Vec3(dist_vx(p.rng), dist_vy(p.rng), dist_vz(p.rng));
                p.geometry.velocity = mesh.transform_direction(
                    p.geometry.position, new_velocity_notransform);
                break;
            }
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
    for (int i = 0; i < ion_processes.size(); i++) {
        p.weight *= std::exp(-dt * ion_processes[i]->rate[cell_id]);
    }
    p.t += dt;

    return event_type;
}

}  // namespace neptune