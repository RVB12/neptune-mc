#include "neptune/finitevolume/advectiondiffusionsolver.hpp"

#include <Eigen/Sparse>

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

    // Construct the Update matrix
    typedef Eigen::SparseMatrix<real_t> SpMat;
    typedef Eigen::Triplet<real_t> Trip;

    std::vector<Trip> tripletList;

    // List of the gamma values in the order of the ghost cells
    std::vector<real_t> gamma_vals;

    index_t n_ghost_cells = 0;

    // Calculate the coefficients for the update matrix
    // for each cell
    for (index_t i = 0; i < n_cells; i++) {
        double diag_coeff = 1.0;

        Vec3 center = mesh.get_center(i);
        real_t volume = mesh.get_volume(i);

        // Loop over the neighbours of the cell
        for (index_t j = 0; j < mesh.get_n_surfaces(i); j++) {
            index_t neighbour_id = mesh.get_neighbour(i, j);

            // note, we assume this normal points inwards
            Vec3 normal = mesh.get_surface_normal(i, j, center);
            // negate the normal to point outwards
            normal = -normal;
            real_t area = mesh.get_surface_area(i, j);

            if (neighbour_id >= 0) {
                // normal cell-cell interaction
                Vec3 neighbour_center = mesh.get_center(neighbour_id);

                real_t neighbour_coeff = 0.0;

                //
                // Diffusion term
                //
                real_t Dr =
                    0.5 * (diffusion_rate[i] + diffusion_rate[neighbour_id]);

                Vec3 dpdx = (neighbour_center - center);
                dpdx /= dpdx.norm2();
                // dpdx *=
                //     (density[neighbour_id] * diffusion_factor[neighbour_id] -
                //      density[i] * diffusion_factor[i]);

                real_t temp =
                    Dr * dot(dpdx, normal) * area * dt / volume;

                diag_coeff += temp * diffusion_factor[i];
                neighbour_coeff -= temp * diffusion_factor[neighbour_id];

                //
                // Advection term (upwind)
                //
                Vec3 advection_vector =
                    0.5 * Vec3(advection_x[i] + advection_x[neighbour_id],
                               advection_y[i] + advection_y[neighbour_id],
                               advection_z[i] + advection_z[neighbour_id]);
                real_t v = dot(advection_vector, normal);
                if (v > 0) {
                    diag_coeff += v * area * dt / volume;
                }
                else {
                    neighbour_coeff += v * area * dt / volume;
                }
                // Set the neighbour coefficient
                tripletList.push_back(Trip(i, neighbour_id, neighbour_coeff));
            }
            else {
                // boundary condition

                // introduce a ghost cell at the cell boundary
                index_t ghost_id = n_cells + n_ghost_cells;
                n_ghost_cells++;
                real_t neighbour_coeff = 0.0;

                //
                // Diffusion term
                //
                Vec3 surface_center = mesh.get_surface_center(i, j);
                real_t Dr = diffusion_rate[i];

                Vec3 dpdx = (surface_center - center);
                real_t dist = dpdx.norm();
                dpdx /= dpdx.norm2();
                // dpdx *=
                //     (density[neighbour_id] * diffusion_factor[neighbour_id] -
                //      density[i] * diffusion_factor[i]);

                real_t temp =
                    Dr * dot(dpdx, normal) * area * dt / volume;

                diag_coeff += temp * diffusion_factor[i];
                neighbour_coeff -= temp * diffusion_factor[i];

                //
                // Advection term (upwind)
                //
                Vec3 advection_vector =
                    Vec3(advection_x[i], advection_y[i], advection_z[i]);
                real_t v = dot(advection_vector, normal);
                if (v > 0) {
                    diag_coeff += v * area * dt / volume;
                }
                else {
                    neighbour_coeff += v * area * dt / volume;
                }
                // Set the neighbour coefficient
                tripletList.push_back(Trip(i, ghost_id, neighbour_coeff));

                // Also add the dynamics for the value of the ghost cell
                // alpha * f - beta * df/dn = gamma
                // => alpha * f - beta * (f - f_ghost) / dist = gamma
                // => (alpha - beta / dist) * f + beta / dist * f_ghost = gamma
                index_t boundary_id = -neighbour_id - 1;
                real_t alpha = boundary_conditions[boundary_id]->alpha;
                real_t beta = boundary_conditions[boundary_id]->beta;
                real_t gamma = boundary_conditions[boundary_id]->gamma;

                tripletList.push_back(Trip(ghost_id, i, (alpha - beta / dist)));
                tripletList.push_back(Trip(ghost_id, ghost_id, beta / dist));
                gamma_vals.push_back(gamma);
            }
        }
        tripletList.push_back(Trip(i, i, diag_coeff));
    }

    SpMat update_matrix(n_cells + n_ghost_cells, n_cells + n_ghost_cells);
    update_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(update_matrix);
    solver.factorize(update_matrix);

    real_t t = 0;
    Eigen::VectorXd rhs(n_cells + n_ghost_cells);
    Eigen::VectorXd solution(n_cells + n_ghost_cells);
    while (t < t_final) {
        // Determine time step
        real_t dt2 = std::min(dt, t_final - t);

        // Construct the right-hand side
        for (index_t i = 0; i < n_cells; i++) {
            rhs[i] = density[i];
        }
        // Add the boundary condition contributions
        for (index_t i = 0; i < n_ghost_cells; i++) {
            rhs[n_cells + i] = gamma_vals[i];
        }

        // Solve the system
        solution = solver.solve(rhs);

        // Update the density
        for (index_t i = 0; i < n_cells; i++) {
            density[i] = solution[i];
        }

        for (index_t i = 0; i < n_cells; i++) {
            density_integral[i] += density[i] * dt2;
        }

        t += dt2;
    }
}

}  // namespace neptune