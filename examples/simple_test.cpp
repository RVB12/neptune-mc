#include <eiron2/geometry/voxelmesh.hpp>
#include <eiron2/simulator/simulator.hpp>
#include <eiron2/common/vector.hpp>
#include <iostream>

using namespace eiron2;

int main()
{
    real_t Lx = 1;
    real_t Ly = 1;
    index_t nx = 100;
    index_t ny = 1;

    VoxelMesh mesh(Lx, Ly, nx, ny);

    index_t n_cells = mesh.n_cells();

    vector<real_t> plasma_dens(n_cells);
    vector<real_t> plasma_temp(n_cells);
    vector<real_t> plasma_vel_x(n_cells);
    vector<real_t> plasma_vel_y(n_cells);

    for (index_t i = 0; i < n_cells; i++) {
        plasma_dens[i] = 1.0;
        plasma_temp[i] = 1.0;
        plasma_vel_x[i] = 0;
        plasma_vel_y[i] = 0;
    }

    PlasmaBackground background(plasma_dens, plasma_temp, plasma_vel_x,
                                plasma_vel_y);

    SimulationParams params;
    index_t particles_per_cell = 100;
    params.n_particles = n_cells * particles_per_cell;
    params.max_steps = 1000;

    vector<real_t> density(n_cells);
    for (index_t i = 0; i < n_cells; i++) {
        density[i] = 0;
    }

    simulate(mesh, background, density, params);

    std::cout << "calculated density" << std::endl;

    for (index_t i = 0; i < n_cells; i++) {
        std::cout << density[i] << std::endl;
    }

    return 0;
}