#include <chrono>
#include <fstream>

#include "neptune/common/json_conversions.hpp"
#include "neptune/simulator/simulator.hpp"

using namespace neptune;

int main(int argc, char** argv)
{
    std::string input_file_name;
    std::string output_file_name;
    if (argc != 3) {
        std::cout << "No input file specified, using default input.json and "
                     "output.json"
                  << std::endl;
        input_file_name = "input.json";
        output_file_name = "output.json";
    }
    else {
        input_file_name = argv[1];
        output_file_name = argv[2];
    }

    std::ifstream input_file(input_file_name);

    if (!input_file.is_open()) {
        std::cerr << "Could not open file " << input_file_name << std::endl;
        return 1;
    }

    std::ofstream output_file(output_file_name);

    if (!output_file.is_open()) {
        std::cerr << "Could not open file " << output_file_name << std::endl;
        return 1;
    }

    nlohmann::json j;

    input_file >> j;

    ExecutionStatistics stats;

    auto start = std::chrono::high_resolution_clock::now();
    std::unique_ptr<Mesh> mesh_ptr = j["mesh"];
    auto end = std::chrono::high_resolution_clock::now();
    stats.grid_building_time =
        std::chrono::duration<double>(end - start).count();

    vector<std::unique_ptr<BGKScatteringProcess>> bgk_processes =
        j["bgk_processes"];

    vector<std::unique_ptr<IonisationProcess>> ion_processes(0);
    if (j.contains("ionisation_processes")) {
        ion_processes =
            j["ionisation_processes"];
    }

    vector<std::unique_ptr<Source>> sources = j["sources"];

    SimulationParams params;
    params.max_steps = j["max_steps"];
    params.seed = j["seed"];

    vector<std::unique_ptr<BoundaryCondition>> boundary_conditions =
        j["boundary_conditions"];

    start = std::chrono::high_resolution_clock::now();
    if (j.contains("t_final")) {
        params.t_final = j["t_final"];
    }
    else {
        params.t_final = -1;
    }

    Estimator estimator(mesh_ptr->n_cells());

    if (j.contains("dt")) {
        real_t dt = j["dt"];
        simulate_kdmc(*mesh_ptr, bgk_processes, ion_processes, sources,
                      estimator, boundary_conditions, params, stats, dt);
    }
    else {
        simulate(*mesh_ptr, bgk_processes, ion_processes, sources, estimator,
                 boundary_conditions, params, stats);
    }

    end = std::chrono::high_resolution_clock::now();

    stats.trace_time = std::chrono::duration<double>(end - start).count();

    stats.total_time = stats.grid_building_time + stats.trace_time;

    nlohmann::json output;

    output["density"] = estimator.density;
    output["momentum_x"] = estimator.momentum_x;
    output["momentum_y"] = estimator.momentum_y;
    output["momentum_z"] = estimator.momentum_z;
    output["energy"] = estimator.energy;

    output["density_final"] = estimator.density_final;
    output["momentum_x_final"] = estimator.momentum_x_final;
    output["momentum_y_final"] = estimator.momentum_y_final;
    output["momentum_z_final"] = estimator.momentum_z_final;
    output["energy_final"] = estimator.energy_final;

    output["stats"] = stats;

    output_file << output.dump(4);

    return 0;
}