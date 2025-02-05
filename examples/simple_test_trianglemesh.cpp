#include <eigrene/common/vector.hpp>
#include <eigrene/geometry/trianglemesh.hpp>
#include <eigrene/physics/linesource.hpp>
#include <eigrene/physics/source.hpp>
#include <eigrene/simulator/simulator.hpp>
#include <iostream>

using namespace eigrene;

int main()
{

    vector<Vec2> vertices(4);
    vertices[0] = Vec2(0, 0);
    vertices[1] = Vec2(1, 0);
    vertices[2] = Vec2(1, 1);
    vertices[3] = Vec2(0, 1);

    vector<TriangleIndices> triangles(2);
    triangles[0] = TriangleIndices{0, 1, 2};
    triangles[1] = TriangleIndices{0, 2, 3};

    // Normally, we would still define the boundaries, but here,
    // we're just going to make everything reflective so it's fine.

    TriangleMesh mesh(vertices, triangles);

    index_t n_cells = mesh.n_cells();

    PlasmaBackground background(n_cells);

    for (index_t i = 0; i < n_cells; i++) {
        background.plasma_dens[i] = 1.0e13;
        background.plasma_temp[i] = 1.0;
        background.plasma_vel_x[i] = 1.0e2;
        background.plasma_vel_y[i] = 0;
    }

    SimulationParams params;
    params.max_steps = 1000;

    vector<real_t> density(n_cells);

    index_t particles_per_cell = 10000;
    LineSource line =
        LineSource(n_cells * particles_per_cell, 1, Vec2(0, 0), Vec2(0, 1));
    vector<std::unique_ptr<Source>> sources(1);
    sources[0] = std::make_unique<LineSource>(line);

    vector<BoundaryCondition> boundary_conditions = {BoundaryCondition::REFLECTIVE};

    ExecutionStatistics stats;
    
    simulate(mesh, background, sources, density, boundary_conditions, params, stats);

    std::cout << "calculated density" << std::endl;

    for (index_t i = 0; i < n_cells; i++) {
        std::cout << density[i] << std::endl;
    }

    return 0;
}