#include <eiron2/geometry/voxelmesh.hpp>
#include <eiron2/simulator/simulator.hpp>
#include <eiron2/common/vector.hpp>
#include <eiron2/physics/source.hpp>
#include <eiron2/physics/linesource.hpp>
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

    PlasmaBackground background(n_cells);

    for (index_t i = 0; i < n_cells; i++) {
        background.plasma_dens[i] = 1.0;
        background.plasma_temp[i] = 1.0;
        background.plasma_vel_x[i] = 0;
        background.plasma_vel_y[i] = 0;
    }

    SimulationParams params;
    params.max_steps = 1000;

    vector<real_t> density(n_cells);
    for (index_t i = 0; i < n_cells; i++) {
        density[i] = 0;
    }

    index_t particles_per_cell = 100;
    LineSource line = LineSource(n_cells * particles_per_cell, Vec2(0, 0), Vec2(0, 1));
    vector<Source*> sources(1);
    sources[0] = &line;

    simulate(mesh, background, sources, density, params);

    std::cout << "calculated density" << std::endl;

    for (index_t i = 0; i < n_cells; i++) {
        std::cout << density[i] << std::endl;
    }

    return 0;
}