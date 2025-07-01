#include "neptune/finitevolume/advectiondiffusionsolver.hpp"

namespace neptune {

void advection_diffusion_solver(
    const VoxelMesh& mesh,
    const vector<std::unique_ptr<FluidBoundaryCondition>>& boundary_conditions,
    const vector<real_t>& advection_x,
    const vector<real_t>& advection_y,
    const vector<real_t>& advection_z,
    const vector<real_t>& diffusion_rate,
    const vector<real_t>& diffusion_factor,
    real_t t_final,
    real_t dt,
    vector<real_t>& density,
    vector<real_t>& density_integral);

void advection_diffusion_solver(
    const TriangleMesh& mesh,
    const vector<std::unique_ptr<FluidBoundaryCondition>>& boundary_conditions,
    const vector<real_t>& advection_x,
    const vector<real_t>& advection_y,
    const vector<real_t>& advection_z,
    const vector<real_t>& diffusion_rate,
    const vector<real_t>& diffusion_factor,
    real_t t_final,
    real_t dt,
    vector<real_t>& density,
    vector<real_t>& density_integral);

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
    // Try to cast the mesh to a VoxelMesh
    const VoxelMesh* voxel_mesh = dynamic_cast<const VoxelMesh*>(&mesh);
    if (voxel_mesh) {
        advection_diffusion_solver(*voxel_mesh, boundary_conditions,
                                   advection_x, advection_y, advection_z,
                                   diffusion_rate, diffusion_factor, t_final,
                                   dt, density, density_integral);
        return;
    }
    // Try to cast the mesh to a TriangleMesh
    const TriangleMesh* triangle_mesh =
        dynamic_cast<const TriangleMesh*>(&mesh);
    if (triangle_mesh) {
        advection_diffusion_solver(*triangle_mesh, boundary_conditions,
                                   advection_x, advection_y, advection_z,
                                   diffusion_rate, diffusion_factor, t_final,
                                   dt, density, density_integral);
        return;
    }
    // If the mesh is neither a VoxelMesh nor a TriangleMesh, throw an error
    throw std::runtime_error(
        "Unsupported mesh type for advection-diffusion solver.");
}

void advection_diffusion_solver(
    const VoxelMesh& mesh,
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

    real_t dx = mesh.cell_size().x;
    real_t dy = mesh.cell_size().y;
    real_t dz = mesh.cell_size().z;

    for (index_t i = 0; i < n_cells; i++) {
        density_integral[i] = 0.0;
    }

    real_t t = 0;
    while (t < t_final) {
        // Determine time step
        real_t dt2 = std::min(dt, t_final - t);

        for (index_t i = 0; i < n_cells; i++) {
            Index3D ind = mesh.linear_to_tuple(i);

            total_flux[i] = 0.0;

            //
            // Calculate the flux on each face
            //
            if (ind.i > 0) {
                real_t area = mesh.get_surface_area(i, 0);
                Index3D ind2(ind.i - 1, ind.j, ind.k);
                index_t i2 = mesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] -
                                  density[i2] * diffusion_factor[i2]) /
                                 (dx)*area;
                // Advection term (upwind)
                real_t v = 0.5 * (advection_x[i] + advection_x[i2]);
                if (v > 0) {
                    total_flux[i] += v * density[i2] * area;
                }
                else {
                    total_flux[i] += v * density[i] * area;
                }
            }
            else {
                // Boundary condition
                real_t area = mesh.get_surface_area(i, 0);
                real_t alpha = boundary_conditions[0]->alpha;
                real_t beta = boundary_conditions[0]->beta;
                real_t gamma = boundary_conditions[0]->gamma;

                real_t ghost =
                    (gamma + density[i] * beta) / (beta - alpha * 0.5 * dx);
                total_flux[i] -= diffusion_factor[i] * diffusion_rate[i] *
                                 (density[i] - ghost) / (0.5 * dx) * area;

                real_t v = advection_x[i];
                if (v > 0) {
                    total_flux[i] -= v * ghost * area;
                }
                else {
                    total_flux[i] -= v * density[i] * area;
                }
            }
            if (ind.i < mesh.nx - 1) {
                real_t area = mesh.get_surface_area(i, 3);
                Index3D ind2(ind.i + 1, ind.j, ind.k);
                index_t i2 = mesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] -
                                  density[i2] * diffusion_factor[i2]) /
                                 (dx)*area;
                // Advection term (upwind)
                real_t v = 0.5 * (advection_x[i] + advection_x[i2]);
                if (v > 0) {
                    total_flux[i] -= v * density[i] * area;
                }
                else {
                    total_flux[i] -= v * density[i2] * area;
                }
            }
            else {
                // Boundary condition
                real_t area = mesh.get_surface_area(i, 3);
                real_t alpha = boundary_conditions[3]->alpha;
                real_t beta = boundary_conditions[3]->beta;
                real_t gamma = boundary_conditions[3]->gamma;

                real_t ghost =
                    (density[i] * beta - gamma) / (beta + alpha * 0.5 * dx);
                total_flux[i] -= diffusion_factor[i] * diffusion_rate[i] *
                                 (density[i] - ghost) / (0.5 * dx) * area;

                real_t v = advection_x[i];
                total_flux[i] -= v * ghost * area;
            }
            if (ind.j > 0) {
                real_t area = mesh.get_surface_area(i, 1);
                Index3D ind2(ind.i, ind.j - 1, ind.k);
                index_t i2 = mesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] -
                                  density[i2] * diffusion_factor[i2]) /
                                 (dy)*area;
                // Advection term (upwind)
                real_t v = 0.5 * (advection_y[i] + advection_y[i2]);
                if (v > 0) {
                    total_flux[i] += v * density[i2] * area;
                }
                else {
                    total_flux[i] += v * density[i] * area;
                }
            }
            else {
                // Boundary condition
                real_t area = mesh.get_surface_area(i, 1);
                real_t alpha = boundary_conditions[1]->alpha;
                real_t beta = boundary_conditions[1]->beta;
                real_t gamma = boundary_conditions[1]->gamma;

                real_t ghost =
                    (gamma + density[i] * beta) / (beta - alpha * 0.5 * dy);
                total_flux[i] -= diffusion_factor[i] * diffusion_rate[i] *
                                 (density[i] - ghost) / (0.5 * dy) * area;

                real_t v = advection_y[i];
                if (v > 0) {
                    total_flux[i] -= v * ghost * area;
                }
                else {
                    total_flux[i] -= v * density[i] * area;
                }
            }
            if (ind.j < mesh.ny - 1) {
                real_t area = mesh.get_surface_area(i, 4);
                Index3D ind2(ind.i, ind.j + 1, ind.k);
                index_t i2 = mesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] -
                                  density[i2] * diffusion_factor[i2]) /
                                 (dy)*area;
                // Advection term (upwind)
                real_t v = 0.5 * (advection_y[i] + advection_y[i2]);
                if (v > 0) {
                    total_flux[i] -= v * density[i] * area;
                }
                else {
                    total_flux[i] -= v * density[i2] * area;
                }
            }
            else {
                // Boundary condition
                real_t area = mesh.get_surface_area(i, 4);
                real_t alpha = boundary_conditions[4]->alpha;
                real_t beta = boundary_conditions[4]->beta;
                real_t gamma = boundary_conditions[4]->gamma;

                real_t ghost =
                    (density[i] * beta - gamma) / (beta + alpha * 0.5 * dy);
                total_flux[i] -= diffusion_factor[i] * diffusion_rate[i] *
                                 (density[i] - ghost) / (0.5 * dy) * area;

                real_t v = advection_y[i];
                total_flux[i] -= v * ghost * area;
            }
            if (ind.k > 0) {
                real_t area = mesh.get_surface_area(i, 2);
                Index3D ind2(ind.i, ind.j, ind.k - 1);
                index_t i2 = mesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] -
                                  density[i2] * diffusion_factor[i2]) /
                                 (dz)*area;
                // Advection term (upwind)
                real_t v = 0.5 * (advection_z[i] + advection_z[i2]);
                if (v > 0) {
                    total_flux[i] += v * density[i2] * area;
                }
                else {
                    total_flux[i] += v * density[i] * area;
                }
            }
            else {
                // Boundary condition
                real_t area = mesh.get_surface_area(i, 2);
                real_t alpha = boundary_conditions[2]->alpha;
                real_t beta = boundary_conditions[2]->beta;
                real_t gamma = boundary_conditions[2]->gamma;

                real_t ghost =
                    (gamma + density[i] * beta) / (beta - alpha * 0.5 * dz);
                total_flux[i] -= diffusion_factor[i] * diffusion_rate[i] *
                                 (density[i] - ghost) / (0.5 * dz) * area;

                real_t v = advection_z[i];
                if (v > 0) {
                    total_flux[i] -= v * ghost * area;
                }
                else {
                    total_flux[i] -= v * density[i] * area;
                }
            }
            if (ind.k < mesh.nz - 1) {
                real_t area = mesh.get_surface_area(i, 5);
                Index3D ind2(ind.i, ind.j, ind.k + 1);
                index_t i2 = mesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] -
                                  density[i2] * diffusion_factor[i2]) /
                                 (dz)*area;
                // Advection term (upwind)
                real_t v = 0.5 * (advection_z[i] + advection_z[i2]);
                if (v > 0) {
                    total_flux[i] -= v * density[i] * area;
                }
                else {
                    total_flux[i] -= v * density[i2] * area;
                }
            }
            else {
                // Boundary condition
                real_t area = mesh.get_surface_area(i, 5);
                real_t alpha = boundary_conditions[5]->alpha;
                real_t beta = boundary_conditions[5]->beta;
                real_t gamma = boundary_conditions[5]->gamma;

                real_t ghost =
                    (density[i] * beta - gamma) / (beta + alpha * 0.5 * dz);
                total_flux[i] -= diffusion_factor[i] * diffusion_rate[i] *
                                 (density[i] - ghost) / (0.5 * dz) * area;

                real_t v = advection_z[i];
                total_flux[i] -= v * ghost * area;
            }
        }

        // Update the density
        for (index_t i = 0; i < n_cells; i++) {
            density[i] = density[i] + dt2 / (dx * dy * dz) * total_flux[i];
        }

        for (index_t i = 0; i < n_cells; i++) {
            density_integral[i] += density[i] * dt2;
        }

        // Update the time
        t += dt;
    }
}

void advection_diffusion_solver(
    const TriangleMesh& mesh,
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
            for (index_t j = 0; j < 3; j++) {
                index_t neighbour_id = mesh.get_neighbour(i, j);

                Vec3 normal = mesh.get_surface_normal(i, j);
                real_t area = mesh.get_surface_area(i, j);

                if (neighbour_id >= 0) {
                    Vec3 neighbour_center = mesh.get_center(neighbour_id);

                    // Diffusion term
                    real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[j]);

                    Vec3 dpdx = (neighbour_center - center);
                    dpdx /= dpdx.norm();
                    dpdx *= (density[i] * diffusion_factor[i] -
                             density[j] * diffusion_factor[j]);

                    total_flux[i] -= Dr * dot(dpdx, normal) * area;

                    // Advection term (upwind)
                    Vec3 advection_vector =
                        0.5 * Vec3(advection_x[i] + advection_x[j],
                                   advection_y[i] + advection_y[j],
                                   advection_z[i] + advection_z[j]);
                    real_t v = dot(advection_vector, normal);
                    if (v > 0) {
                        total_flux[i] += v * density[j] * area;
                    }
                    else {
                        total_flux[i] += v * density[i] * area;
                    }
                }
                else {
                    // TODO: handle boundary conditions
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
        t += dt;
    }
}

}  // namespace neptune