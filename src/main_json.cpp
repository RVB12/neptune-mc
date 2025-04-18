#include <chrono>
#include <fstream>

#include "eigrene/common/json_conversions.hpp"
#include "eigrene/simulator/simulator.hpp"

using namespace eigrene;

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

    PlasmaBackground background = j["background"];

    vector<std::unique_ptr<Source>> sources = j["sources"];

    SimulationParams params;
    params.max_steps = j["max_steps"];
    params.seed = j["seed"];

    vector<BoundaryCondition> boundary_conditions = j["boundary_conditions"];

    start = std::chrono::high_resolution_clock::now();
    if (j.contains("t_final")){
        params.t_final = j["t_final"];
    }
    else {
        params.t_final = -1;
    }

    Estimator estimator(mesh_ptr->n_cells());

    if( j.contains("dt")){
        real_t dt = j["dt"];
        simulate_kdmc2(*mesh_ptr, background, sources, estimator,
                      boundary_conditions, params, stats, dt);
    } else {
        simulate(*mesh_ptr, background, sources, estimator, boundary_conditions,
                 params, stats);
    }

    end = std::chrono::high_resolution_clock::now();

    stats.trace_time = std::chrono::duration<double>(end - start).count();

    stats.total_time = stats.grid_building_time + stats.trace_time;

    nlohmann::json output;

    output["density"] = estimator.density;

    output["stats"] = stats;

    output_file << output.dump(4);

    return 0;
}