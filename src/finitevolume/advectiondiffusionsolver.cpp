#include "neptune/finitevolume/advectiondiffusionsolver.hpp"

#include <Eigen/Sparse>

namespace neptune {

std::pair<Eigen::SparseMatrix<real_t>, std::vector<real_t>>
construct_update_matrix(
    const Mesh& mesh,
    const vector<std::unique_ptr<FluidBoundaryCondition>>& boundary_conditions,
    const vector<real_t>& advection_x,
    const vector<real_t>& advection_y,
    const vector<real_t>& advection_z,
    const vector<real_t>& diffusion_rate,
    const vector<real_t>& diffusion_factor,
    const vector<real_t>& ionisation_rate,
    real_t dt)
{
    const index_t n_cells = mesh.n_cells();

    typedef Eigen::Triplet<real_t> Trip;

    std::vector<Trip> tripletList;

    std::vector<real_t> gamma_vals;

    index_t n_ghost_cells = 0;

    // Calculate the coefficients for the update matrix
    // for each cell
    for (index_t i = 0; i < n_cells; i++) {
        double diag_coeff = 1.0;

        Vec3 center = mesh.get_center(i);
        real_t volume = mesh.get_volume(i);

        // Add the ionisation term
        diag_coeff += ionisation_rate[i] * dt;

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

                real_t temp = Dr * dot(dpdx, normal) * area * dt / volume;

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
                index_t boundary_id = -neighbour_id - 1;
                FluidBoundaryCondition* bc =
                    boundary_conditions[-neighbour_id - 1].get();

                // Check which type of boundary condition we have
                if (dynamic_cast<NoFluxFluidBoundaryCondition*>(bc)) {
                    // No flux boundary condition, we can just ignore the
                    // neighbour and not add any terms to the diagonal
                }
                else {
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
                    //     (density[neighbour_id] *
                    //     diffusion_factor[neighbour_id] -
                    //      density[i] * diffusion_factor[i]);

                    real_t temp = Dr * dot(dpdx, normal) * area * dt / volume;

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

                    // Now we also need to add some equations to determine the
                    // value of the ghost cell
                    if (dynamic_cast<AbsorbingFluidBoundaryCondition*>(bc)) {
                        // Absorbing boundary condition, we can just set the
                        // ghost cell value to zero
                        tripletList.push_back(Trip(ghost_id, ghost_id, 1.0));
                        gamma_vals.push_back(0.0);
                        continue;
                    }
                    else if (dynamic_cast<RobinFluidBoundaryCondition*>(bc)) {
                        // Robin boundary condition, we will add the equation
                        real_t alpha =
                            dynamic_cast<RobinFluidBoundaryCondition*>(bc)
                                ->alpha;
                        real_t beta =
                            dynamic_cast<RobinFluidBoundaryCondition*>(bc)
                                ->beta;
                        real_t gamma =
                            dynamic_cast<RobinFluidBoundaryCondition*>(bc)
                                ->gamma;

                        if (beta == 0) {
                            tripletList.push_back(
                                Trip(ghost_id, ghost_id, alpha));
                            gamma_vals.push_back(gamma);
                        }
                        else {
                            tripletList.push_back(
                                Trip(ghost_id, i, (alpha - beta / dist)));
                            tripletList.push_back(
                                Trip(ghost_id, ghost_id, beta / dist));
                            gamma_vals.push_back(gamma);
                        }
                    }
                }
            }
        }
        tripletList.push_back(Trip(i, i, diag_coeff));
    }

    Eigen::SparseMatrix<real_t> update_matrix(n_cells + n_ghost_cells,
                                              n_cells + n_ghost_cells);
    update_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    return {update_matrix, gamma_vals};
}

void advection_diffusion_solver(
    const Mesh& mesh,
    const vector<std::unique_ptr<FluidBoundaryCondition>>& boundary_conditions,
    const vector<real_t>& advection_x,
    const vector<real_t>& advection_y,
    const vector<real_t>& advection_z,
    const vector<real_t>& diffusion_rate,
    const vector<real_t>& diffusion_factor,
    const vector<real_t>& ionisation_rate,
    real_t t_final,
    real_t dt,
    vector<real_t>& density,
    vector<real_t>& density_integral,
    vector<real_t>& momentum_x_integral,
    vector<real_t>& momentum_y_integral,
    vector<real_t>& momentum_z_integral,
    vector<real_t>& energy_integral,
    vector<real_t>& energy_x_integral)
{
    index_t n_cells = mesh.n_cells();
    vector<real_t> total_flux(n_cells, 0.0);

    for (index_t i = 0; i < n_cells; i++) {
        density_integral[i] = 0.0;
        momentum_x_integral[i] = 0.0;
        momentum_y_integral[i] = 0.0;
        momentum_z_integral[i] = 0.0;
        energy_integral[i] = 0.0;
        energy_x_integral[i] = 0.0;
    }

    auto [update_matrix, gamma_vals] = construct_update_matrix(
        mesh, boundary_conditions, advection_x, advection_y, advection_z,
        diffusion_rate, diffusion_factor, ionisation_rate, dt);

    size_t n_ghost_cells = gamma_vals.size();

    Eigen::SparseLU<Eigen::SparseMatrix<real_t>> solver;
    solver.analyzePattern(update_matrix);
    solver.factorize(update_matrix);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to factorize the update matrix");
    }

    real_t t = 0;
    Eigen::VectorXd rhs(n_cells + n_ghost_cells);
    Eigen::VectorXd solution(n_cells + n_ghost_cells);
    while (t < t_final) {
        // Determine time step and potentially recalculate the update matrix
        real_t dt2;
        if (t + dt > t_final) {
            dt2 = t_final - t;
            auto [new_update_matrix, new_gamma_vals] = construct_update_matrix(
                mesh, boundary_conditions, advection_x, advection_y,
                advection_z, diffusion_rate, diffusion_factor, ionisation_rate,
                dt2);
            update_matrix = new_update_matrix;
            gamma_vals = new_gamma_vals;
            solver.factorize(update_matrix);

            if (solver.info() != Eigen::Success) {
                throw std::runtime_error(
                    "Failed to factorize the update matrix");
            }
        }
        else {
            dt2 = dt;
        }

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

        // Update the integrals
        for (index_t i = 0; i < n_cells; i++) {
            density_integral[i] += density[i] * dt2;

            // TODO: this differential form is not great, we should do this in a
            // similar flux style
            if (i < n_cells - 1) {
                real_t dx = mesh.get_center(i + 1).x - mesh.get_center(i).x;
                real_t dy = mesh.get_center(i + 1).y - mesh.get_center(i).y;
                real_t dz = mesh.get_center(i + 1).z - mesh.get_center(i).z;
                momentum_x_integral[i] +=
                    (density[i] * advection_x[i] -
                     diffusion_rate[i] *
                         (density[i + 1] * diffusion_factor[i + 1] -
                          density[i] * diffusion_factor[i]) /
                         dx) *
                    dt2;
                momentum_y_integral[i] +=
                    (density[i] * advection_y[i] -
                     diffusion_rate[i] *
                         (density[i + 1] * diffusion_factor[i + 1] -
                          density[i] * diffusion_factor[i]) /
                         dy) *
                    dt2;
                momentum_z_integral[i] +=
                    (density[i] * advection_z[i] -
                     diffusion_rate[i] *
                         (density[i + 1] * diffusion_factor[i + 1] -
                          density[i] * diffusion_factor[i]) /
                         dz) *
                    dt2;

                energy_x_integral[i] +=
                    (0.5 * density[i] *
                         (advection_x[i] * advection_x[i] +
                          diffusion_factor[i]) -
                     diffusion_rate[i] *
                         (density[i + 1] * advection_x[i + 1] *
                              diffusion_factor[i + 1] -
                          density[i] * advection_x[i] * diffusion_factor[i]) /
                         dx) *
                    dt2;
                energy_integral[i] +=
                    (0.5 * density[i] *
                         (advection_x[i] * advection_x[i] +
                          advection_y[i] * advection_y[i] +
                          advection_z[i] * advection_z[i] +
                          3.0 * diffusion_factor[i]) -
                     diffusion_rate[i] *
                         (density[i + 1] * advection_x[i + 1] *
                              diffusion_factor[i + 1] -
                          density[i] * advection_x[i] * diffusion_factor[i]) /
                         dx) *
                    dt2;
            }
            else {
                real_t dx = mesh.get_center(i - 1).x - mesh.get_center(i).x;
                real_t dy = mesh.get_center(i - 1).y - mesh.get_center(i).y;
                real_t dz = mesh.get_center(i - 1).z - mesh.get_center(i).z;
                momentum_x_integral[i] +=
                    (density[i] * advection_x[i] -
                     diffusion_rate[i] *
                         (density[i - 1] * diffusion_factor[i - 1] -
                          density[i] * diffusion_factor[i]) /
                         dx) *
                    dt2;
                momentum_y_integral[i] +=
                    (density[i] * advection_y[i] -
                     diffusion_rate[i] *
                         (density[i - 1] * diffusion_factor[i - 1] -
                          density[i] * diffusion_factor[i]) /
                         dy) *
                    dt2;
                momentum_z_integral[i] +=
                    (density[i] * advection_z[i] -
                     diffusion_rate[i] *
                         (density[i - 1] * diffusion_factor[i - 1] -
                          density[i] * diffusion_factor[i]) /
                         dz) *
                    dt2;

                energy_x_integral[i] +=
                    (0.5 * density[i] *
                         (advection_x[i] * advection_x[i] +
                          diffusion_factor[i]) -
                     diffusion_rate[i] *
                         (density[i - 1] * advection_x[i - 1] *
                              diffusion_factor[i - 1] -
                          density[i] * advection_x[i] * diffusion_factor[i]) /
                         dx) *
                    dt2;
                energy_integral[i] +=
                    (0.5 * density[i] *
                         (advection_x[i] * advection_x[i] +
                          advection_y[i] * advection_y[i] +
                          advection_z[i] * advection_z[i] +
                          3.0 * diffusion_factor[i]) -
                     diffusion_rate[i] *
                         (density[i - 1] * advection_x[i - 1] *
                              diffusion_factor[i - 1] -
                          density[i] * advection_x[i] * diffusion_factor[i]) /
                         dx) *
                    dt2;
            }
        }

        t += dt2;
    }

    // Multiple energy integrals by 2 * constants::m_p / constants::eV_to_Joules
    // to change units
    for (index_t i = 0; i < n_cells; i++) {
        energy_x_integral[i] *= 2.0 * constants::m_p / constants::eV_to_Joules;
        energy_integral[i] *= 2.0 * constants::m_p / constants::eV_to_Joules;
    }
}

}  // namespace neptune