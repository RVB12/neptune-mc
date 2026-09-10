#include "neptune/physics/update_particle.hpp"

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/geometry/cart_to_tor.hpp"
#include "neptune/physics/background.hpp"
#include "neptune/physics/physicalconstants.hpp"
namespace neptune {

EventType update_particle(
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ion_collisions,
    const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions,
    Particle& p,
    real_t t_final)
{
    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);
    std::exponential_distribution<real_t> r_exp(1.0);

    index_t cell_id = p.geometry.cell_id;  // Current cell id of the particle

    // Calculate total collision rate (charge exchange + neutral-neutral)

    real_t collision_rate = 0;
    for (int i = 0; i < cx_collisions.size(); i++) {
        collision_rate += cx_collisions[i]->collision_rate(p, background);
    }
    // time till collision
    real_t tc = p.u / collision_rate;

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
            for (int i = 0; i < ion_collisions.size(); i++) {
                // TODO: only consider collisions for the current particle
                // species
                real_t R_i = ion_collisions[i]->collision_rate(p, background);
                p.weight *= std::exp(-time_left * R_i);
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

        if (cx_collisions.size() == 1) {
            Vec3 new_velocity_notransform =
                cx_collisions[0]->sample_post_collision_velocity(p, background);
            if (mesh.is_toroidal()) {
                // If the mesh uses toroidal coordinates, we need to
                // transform the velocity
                new_velocity_notransform = transform_tor_dir_to_cart(
                    p.geometry.position, new_velocity_notransform);
            }
            else {
                p.geometry.velocity = new_velocity_notransform;
            }
            p.u = r_exp(p.rng);
        }
        else {
            if (cx_collisions.size() == 1) {
                // Only one collision type, we can directly sample the
                // post-collision velocity
                Vec3 new_velocity_notransform =
                    cx_collisions[0]->sample_post_collision_velocity(
                        p, background);
                if (mesh.is_toroidal()) {
                    // If the mesh uses toroidal coordinates, we need to
                    // transform the velocity
                    new_velocity_notransform = transform_tor_dir_to_cart(
                        p.geometry.position, new_velocity_notransform);
                }
                else {
                    p.geometry.velocity = new_velocity_notransform;
                }
                p.u = r_exp(p.rng);
            }
            else {
                // A collision has occurred, but we don't know which one
                // draw a new random number to determine that
                real_t eps = r_unif(p.rng);
                real_t cumulative_rate = 0;
                for (int i = 0; i < cx_collisions.size(); i++) {
                    cumulative_rate +=
                        cx_collisions[i]->collision_rate(p, background);
                    if (eps <= cumulative_rate / collision_rate) {
                        // We have a collision with the background
                        // Generate new velocity from plasma background (BGK
                        // operator)
                        Vec3 new_velocity_notransform =
                            cx_collisions[i]->sample_post_collision_velocity(
                                p, background);
                        if (mesh.is_toroidal()) {
                            // If the mesh uses toroidal coordinates, we need to
                            // transform the velocity
                            new_velocity_notransform =
                                transform_tor_dir_to_cart(
                                    p.geometry.position,
                                    new_velocity_notransform);
                        }
                        else {
                            p.geometry.velocity = new_velocity_notransform;
                        }
                        p.u = r_exp(p.rng);
                        break;
                    }
                }
            }
        }

        // Reset boundary id (we are no longer at a boundary)
        p.geometry.boundary_id = -1;

        dt = tc;
    }
    else {
        event_type = EventType::CELL_CROSSING;
        dt = crossing.t;
        p.u -=
            collision_rate * dt;  // Update the scaled time till next collision
        if (crossing.boundary_condition_id == -1) {
            // Regular grid crossing, move the particle to the next
            // cell
            p.geometry.position = crossing.position;
            p.geometry.cell_id = crossing.next_cell_id;
            p.geometry.boundary_id = crossing.next_surface_id;
        }
        else {
            // There is a boundary condition to be applied

            // First, move the particle to the boundary
            p.geometry.position = crossing.position;
            p.geometry.boundary_id = crossing.surface_id;

            boundary_conditions[crossing.boundary_condition_id]->apply(p, mesh, crossing);
        }
    }

    // Ionisation
    for (int i = 0; i < ion_collisions.size(); i++) {
        p.weight *=
            std::exp(-dt * ion_collisions[i]->collision_rate(p, background));
    }
    p.t += dt;

    return event_type;
}

}  // namespace neptune