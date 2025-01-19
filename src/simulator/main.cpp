#include "eiron2/simulator/simulator.hpp"
#include "eiron2/common/json_conversions.hpp"

#include <fstream>

using namespace eiron2;

int main(int argc, char** argv){

    if(argc != 2){
        std::cerr << "Usage: " << argv[0] << " <input_file.json>" << std::endl;
        return 1;
    }

    std::ifstream input_file(argv[1]);

    if(!input_file.is_open()){
        std::cerr << "Could not open file " << argv[1] << std::endl;
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

    simulate(*mesh_ptr, background, sources, density, params);

    std::cout << "calculated density" << std::endl;

    for (index_t i = 0; i < density.size(); i++) {
        std::cout << density[i] << std::endl;
    }

    return 0;
}