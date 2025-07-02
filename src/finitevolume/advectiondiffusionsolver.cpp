#include "neptune/finitevolume/advectiondiffusionsolver.hpp"

namespace neptune {

void advection_diffusion_solver(
    const Mesh& mesh,
    const vector<std::unique_ptr<FluidBoundaryCondition>>& boundary_conditions,
    const vector<real_t>& advection_x,
    const vector<real_t>& advection_y,
    const vector<real_t>& advection_z,
    const vector<real_t>& diffusion_rate,
    const vector<real_t>& diffusion_factor,
    real_t t_final,
    real_t dt,
    vector<real_t>& density,
    vector<real_t>& density_integral)
{
    index_t n_cells = mesh.n_cells();
    vector<real_t> total_flux(n_cells, 0.0);

    for (index_t i = 0; i < n_cells; i++) {
        density_integral[i] = 0.0;
    }

    real_t t = 0;
    while (t < t_final) {
        // Determine time step
        real_t dt2 = std::min(dt, t_final - t);

        // Calculate the flux on each face
        for (index_t i = 0; i < n_cells; i++) {
            total_flux[i] = 0.0;

            Vec3 center = mesh.get_center(i);

            // Loop over the neighbours of the cell
            for (index_t j = 0; j < mesh.get_n_surfaces(i); j++) {
                index_t neighbour_id = mesh.get_neighbour(i, j);

                // note, we assume this normal points inwards
                Vec3 normal = mesh.get_surface_normal(i, j);
                // negate the normal to point outwards
                normal = -normal;
                real_t area = mesh.get_surface_area(i, j);

                if (neighbour_id >= 0) {
                    Vec3 neighbour_center = mesh.get_center(neighbour_id);

                    // Diffusion term
                    real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[neighbour_id]);

                    Vec3 dpdx = (neighbour_center - center);
                    dpdx /= dpdx.norm2();
                    dpdx *= (density[neighbour_id] * diffusion_factor[neighbour_id] -
                             density[i] * diffusion_factor[i]);

                    total_flux[i] += Dr * dot(dpdx, normal) * area;

                    // Advection term (upwind)
                    Vec3 advection_vector =
                        0.5 * Vec3(advection_x[i] + advection_x[neighbour_id],
                                   advection_y[i] + advection_y[neighbour_id],
                                   advection_z[i] + advection_z[neighbour_id]);
                    real_t v = dot(advection_vector, normal);
                    if (v > 0) {
                        total_flux[i] -= v * density[i] * area;
                    }
                    else {
                        total_flux[i] -= v * density[neighbour_id] * area;
                    }
                }
                else {
                    index_t boundary_id = -neighbour_id - 1;
                    real_t alpha = boundary_conditions[boundary_id]->alpha;
                    real_t beta = boundary_conditions[boundary_id]->beta;
                    real_t gamma = boundary_conditions[boundary_id]->gamma;


                    Vec3 surface_center = mesh.get_surface_center(i,j);

                    Vec3 dir = surface_center - center;
                    real_t dist = dir.norm();

                    real_t ghost = (dist * gamma + density[i] * beta) /
                                   (beta + alpha * dist);
                    real_t dpdn = (ghost - density[i]) / dist;
                    total_flux[i] += diffusion_factor[i] * diffusion_rate[i] *
                                     dpdn * area;

                    Vec3 advection_vector = Vec3(advection_x[i], advection_y[i], advection_z[i]);
                    real_t v = dot(advection_vector, normal);
                    total_flux[i] -= v * ghost * area;
                }
            }
        }

        // Update the density
        for (index_t i = 0; i < n_cells; i++) {
            density[i] = density[i] + dt2 / mesh.get_volume(i) * total_flux[i];
        }

        for (index_t i = 0; i < n_cells; i++) {
            density_integral[i] += density[i] * dt2;
        }

        // Update the time
        t += dt2;
    }
}

}  // namespace neptune