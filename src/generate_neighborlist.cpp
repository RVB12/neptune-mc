#include <chrono>
#include <fstream>

#include "neptune/common/json_conversions.hpp"
#include "neptune/common/version.hpp"
#include "neptune/simulator/simulator.hpp"
#include "neptune/geometry/trianglemesh.hpp"

using namespace neptune;

int main(int argc, char** argv)
{
    //
    // Parse command line and read input file
    //
    std::string input_file_name;
    if (argc >= 2) {
        std::string argument = argv[1];
        if (argument == "-v" || argument == "--version") {
            std::cout << "Neptune version: " << NEPTUNE_VERSION << std::endl;
            return 0;
        }
    } else {
        std::cout << "No input file specified, using default input.h5"
                  << std::endl;
        input_file_name = "input.h5";
    }

    std::ifstream input_file(input_file_name);

    if (!input_file.is_open()) {
        std::cerr << "Could not open file " << input_file_name << std::endl;
        return 1;
    }

    //
    // Parse input file
    //
    nlohmann::json j;

    input_file >> j;

    auto start = std::chrono::high_resolution_clock::now();
    std::unique_ptr<Mesh> mesh_ptr = j["mesh"];
    auto end = std::chrono::high_resolution_clock::now();

    TriangleMesh* triangle_mesh =
        dynamic_cast<TriangleMesh*>(mesh_ptr.get());

    if (!triangle_mesh) {
        std::cerr << "Mesh is not a TriangleMesh" << std::endl;
        return 1;
    }

    //
    // Write output to file
    //
    // nlohmann::json output;
    
    // output["neighbor_list"] = triangle_mesh->get_neighbor_list();

    // output_file << output.dump(4);

    return 0;
}