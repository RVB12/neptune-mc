#include "eiron2/physics/update_particle.hpp"

namespace eiron2 {

void update_particle(const Mesh& mesh, const PlasmaBackground& background, Particle& p)
{

    // TODO: Move these constants to a more appropriate place
    real_t T_i = 40;
    real_t T_e = T_i;
    // real_t T_n = T_i;
    real_t K_cx = 3.2e-15 * std::sqrt(T_i / 0.026);
    real_t K_i = 2.0e-13 / (6.0 + (T_e / 13.6)) * std::sqrt(T_e / 13.6) *
                 std::exp(-13.6 / T_e);

    std::uniform_real_distribution<real_t> r_unif(0, 1);
    std::normal_distribution<real_t> r_normal(0, 1);

    index_t cell_id =
        p.geometry.cell_id;  // Current cell id of the particle

    // random number to determine time till collision
    real_t eps = r_unif(p.rng);
    // time till collision
    real_t tc = -std::log(eps) / (K_cx * background.plasma_dens[cell_id]);

    // Calculate the next cell crossing
    CellCrossing crossing = mesh.get_next_crossing(p.geometry);

    // Check which one is earlier, charge exchange or cell crossing
    real_t dt;
    if (tc < crossing.t) {
        // Move the particle until the collision
        p.geometry.position += tc * p.geometry.velocity;
        // Generate new velocity from plasma background (BGK
        // operator)
        real_t stddev = std::sqrt(background.plasma_temp[cell_id]);
        std::normal_distribution<real_t> dist_vx(
            background.plasma_vel_x[cell_id], stddev);
        std::normal_distribution<real_t> dist_vy(
            background.plasma_vel_y[cell_id], stddev);
        p.geometry.velocity =
            Vec2(dist_vx(p.rng), dist_vy(p.rng));

        dt = tc;
    }
    else {
        dt = crossing.t;
        if (crossing.next_cell_id != -1) {
            // Regular grid crossing, move the particle to the next
            // cell
            p.geometry.position = crossing.position;
            p.geometry.cell_id = crossing.next_cell_id;
        }
        else {
            // Boundary crossing
            // For now, we just reflect particles at boundaries
            // 0,1,2,3,4 and absorb particles at boundary 3 (x ==
            // Xl)

            // Todo: Implement proper boundary conditions
            // maybe implement surface norm in the mesh

            p.geometry.position = crossing.position;
            if (crossing.boundary_id == 0) {
                p.geometry.velocity.x = -p.geometry.velocity.x;
            }
            else if (crossing.boundary_id == 1 || crossing.boundary_id == 3) {
                p.geometry.velocity.y = -p.geometry.velocity.y;
            }
            else if (crossing.boundary_id == 2) {
                p.weight = 0;
            }
            else {
                throw std::runtime_error("Unknown boundary id");
            }
        }
    }

    // Ionisation
    p.weight =
        p.weight * std::exp(-dt * K_i * background.plasma_dens[cell_id]);

    p.t += dt;
}

}  // namespace eiron2