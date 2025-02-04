#include "eiron2/physics/update_particle.hpp"

#include "eiron2/physics/physicalconstants.hpp"

namespace eiron2 {

void update_particle(const Mesh& mesh,
                     const PlasmaBackground& background,
                     const vector<BoundaryCondition>& boundary_conditions,
                     Particle& p,
                     index_t& num_crossings,
                     index_t& num_collisions)
{
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

    // Check which one is earlier, charge exchange or cell crossing
    real_t dt;
    if (tc < crossing.t) {
        num_collisions++;
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
        p.geometry.velocity = Vec2(dist_vx(p.rng), dist_vy(p.rng));

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
                Vec2 normal =
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

    // Ionisation
    p.weight =
        p.weight *
        std::exp(-dt *
                 constants::calculate_K_i(background.plasma_temp[cell_id]) *
                 background.plasma_dens[cell_id]);

    p.t += dt;
}

}  // namespace eiron2