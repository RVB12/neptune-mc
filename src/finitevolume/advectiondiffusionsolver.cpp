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

    // Cast to voxel mesh
    const VoxelMesh& vmesh = dynamic_cast<const VoxelMesh&>(mesh);
    real_t dx = vmesh.cell_size().x;
    real_t dy = vmesh.cell_size().y;
    real_t dz = vmesh.cell_size().z;

    for (index_t i = 0; i < n_cells; i++) {
        density_integral[i] = 0.0;
    }

    real_t t = 0;
    while (t < t_final) {
        // Determine time step
        real_t dt2 = std::min(dt, t_final - t);

        for (index_t i = 0; i < n_cells; i++) {
            Index3D ind = vmesh.linear_to_tuple(i);

            total_flux[i] = 0.0;

            //
            // Calculate the flux on each face
            //
            if (ind.i > 0) {
                Index3D ind2(ind.i - 1, ind.j, ind.k);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] - density[i2] * diffusion_factor[i2]) / (dx);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_x[i] + advection_x[i2]);
                if (v > 0) {
                    total_flux[i] += v * density[i2];
                }
                else {
                    total_flux[i] += v * density[i];
                }
            }
            else {
                // Boundary condition
                real_t alpha = boundary_conditions[0]->alpha;
                real_t beta = boundary_conditions[0]->beta;
                real_t gamma = boundary_conditions[0]->gamma;

                // TODO: double check this
                total_flux[i] += (gamma - alpha * density[i]) /
                                 (beta + alpha * dx / (2 * diffusion_rate[i]*diffusion_factor[i])) * dy *
                                 dz;
            }
            if (ind.i < vmesh.nx - 1) {
                Index3D ind2(ind.i + 1, ind.j, ind.k);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] - density[i2] * diffusion_factor[i2]) / (dx);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_x[i] + advection_x[i2]);
                if (v > 0) {
                    total_flux[i] -= v * density[i];
                }
                else {
                    total_flux[i] -= v * density[i2];
                }
            }
            else {
                // Boundary condition
                real_t alpha = boundary_conditions[3]->alpha;
                real_t beta = boundary_conditions[3]->beta;
                real_t gamma = boundary_conditions[3]->gamma;

                // TODO: double check this
                total_flux[i] += (gamma - alpha * density[i]) /
                                 (beta + alpha * dx / (2 * diffusion_rate[i]*diffusion_factor[i])) * dy *
                                 dz;
            }
            if (ind.j > 0) {
                Index3D ind2(ind.i, ind.j - 1, ind.k);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] - density[i2] * diffusion_factor[i2]) / (dy);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_y[i] + advection_y[i2]);
                if (v > 0) {
                    total_flux[i] += v * density[i2];
                }
                else {
                    total_flux[i] += v * density[i];
                }
            }
            if (ind.j < vmesh.ny - 1) {
                Index3D ind2(ind.i, ind.j + 1, ind.k);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] - density[i2] * diffusion_factor[i2]) / (dy);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_y[i] + advection_y[i2]);
                if (v > 0) {
                    total_flux[i] -= v * density[i];
                }
                else {
                    total_flux[i] -= v * density[i2];
                }
            }
            if (ind.k > 0) {
                Index3D ind2(ind.i, ind.j, ind.k - 1);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] - density[i2] * diffusion_factor[i2]) / (dz);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_z[i] + advection_z[i2]);
                if (v > 0) {
                    total_flux[i] += v * density[i2];
                }
                else {
                    total_flux[i] += v * density[i];
                }
            }
            if (ind.k < vmesh.nz - 1) {
                Index3D ind2(ind.i, ind.j, ind.k + 1);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                real_t Dr = 0.5 * (diffusion_rate[i] + diffusion_rate[i2]);
                total_flux[i] -= Dr *
                                 (density[i] * diffusion_factor[i] - density[i2] * diffusion_factor[i2]) / (dz);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_z[i] + advection_z[i2]);
                if (v > 0) {
                    total_flux[i] -= v * density[i];
                }
                else {
                    total_flux[i] -= v * density[i2];
                }
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

}  // namespace neptune