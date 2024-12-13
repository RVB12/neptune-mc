#include <eiron2/geometry/voxelmesh.hpp>
#include <eiron2/simulator/simulator.hpp>
#include <iostream>
#include <vector>

using namespace eiron2;

int main()
{
    double Lx = 1;
    double Ly = 1;
    double Lz = 1;
    int nx = 20;
    int ny = 1;
    int nz = 1;

    VoxelMesh mesh(Lx, Ly, Lz, nx, ny, nz);

    int n_cells = mesh.n_cells();

    std::vector<double> plasma_dens(n_cells);
    std::vector<double> plasma_temp(n_cells);
    std::vector<double> plasma_vel_x(n_cells);
    std::vector<double> plasma_vel_y(n_cells);
    std::vector<double> plasma_vel_z(n_cells);

    for (int i = 0; i < n_cells; i++) {
        plasma_dens[i] = 1.0;
        plasma_temp[i] = 1.0;
        plasma_vel_x[i] = 0;
        plasma_vel_y[i] = 0;
        plasma_vel_z[i] = 0;
    }

    PlasmaBackground background{plasma_dens, plasma_temp, plasma_vel_x,
                                plasma_vel_y, plasma_vel_z};

    SimulationParams params;
    params.n_particles = 10000;
    params.max_steps = 10000;

    std::vector<double> density(n_cells);
    for (int i = 0; i < n_cells; i++) {
        density[i] = 0;
    }

    simulate(mesh, background, density, params);

    std::cout << "calculated density" << std::endl;

    for (int i = 0; i < n_cells; i++) {
        std::cout << density[i] << std::endl;
    }

    return 0;
}