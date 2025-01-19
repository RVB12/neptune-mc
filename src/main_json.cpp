#include <fstream>

#include "eiron2/common/json_conversions.hpp"
#include "eiron2/simulator/simulator.hpp"

using namespace eiron2;

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

    std::unique_ptr<Mesh> mesh_ptr = j["mesh"];

    PlasmaBackground background = j["background"];

    vector<std::unique_ptr<Source>> sources = j["sources"];

    SimulationParams params;
    params.max_steps = j["max_steps"];
    params.seed = j["seed"];

    vector<real_t> density(mesh_ptr->n_cells());

    vector<BoundaryCondition> boundary_conditions = j["boundary_conditions"];

    simulate(*mesh_ptr, background, sources, density, boundary_conditions,
             params);

    nlohmann::json output;

    output["density"] = density;

    output_file << output.dump(4);

    return 0;
}