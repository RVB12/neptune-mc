#include "neptune/finitevolume/advectiondiffusionsolver.hpp"

namespace neptune {

void advection_diffusion_solver(const Mesh& mesh,
                                const vector<real_t>& advection_x,
                                const vector<real_t>& advection_y,
                                const vector<real_t>& advection_z,
                                const vector<real_t>& diffusion,
                                real_t t_final,
                                real_t dt,
                                vector<real_t>& density)
{
    index_t n_cells = mesh.n_cells();
    vector<real_t> total_flux(n_cells, 0.0);

    // Cast to voxel mesh
    const VoxelMesh& vmesh = dynamic_cast<const VoxelMesh&>(mesh);
    real_t dx = vmesh.cell_size().x;
    real_t dy = vmesh.cell_size().y;
    real_t dz = vmesh.cell_size().z;

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
                total_flux[i] -= 0.5 * (diffusion[i] + diffusion[i2]) *
                                  (density[i] - density[i2]) / (dx);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_x[i] + advection_x[i2]);
                if (v > 0) {
                    total_flux[i] += v * density[i2];
                }
                else {
                    total_flux[i] += v * density[i];
                }
            }
            if (ind.i < vmesh.nx - 1) {
                Index3D ind2(ind.i + 1, ind.j, ind.k);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                total_flux[i] -= 0.5 * (diffusion[i] + diffusion[i2]) *
                                  (density[i] - density[i2]) / (dx);
                // Advection term (upwind)
                real_t v = 0.5 * (advection_x[i] + advection_x[i2]);
                if (v > 0) {
                    total_flux[i] -= v * density[i];
                }
                else {
                    total_flux[i] -= v * density[i2];
                }
            }
            if (ind.j > 0) {
                Index3D ind2(ind.i, ind.j - 1, ind.k);
                index_t i2 = vmesh.tuple_to_linear(ind2);
                // Diffusion term
                total_flux[i] -= 0.5 * (diffusion[i] + diffusion[i2]) *
                                  (density[i] - density[i2]) / (dy);
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
                total_flux[i] -= 0.5 * (diffusion[i] + diffusion[i2]) *
                                  (density[i] - density[i2]) / (dy);
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
                total_flux[i] -= 0.5 * (diffusion[i] + diffusion[i2]) *
                                  (density[i] - density[i2]) / (dz);
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
                total_flux[i] -= 0.5 * (diffusion[i] + diffusion[i2]) *
                                  (density[i] - density[i2]) / (dz);
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

        // Update the time
        t += dt;
    }
}

}  // namespace neptune