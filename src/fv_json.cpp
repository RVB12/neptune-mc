#include <chrono>
#include <fstream>

#include "neptune/common/json_conversions.hpp"
#include "neptune/finitevolume/advectiondiffusionsolver.hpp"

using namespace neptune;

int main(int argc, char** argv)
{
    //
    // Parse command line and read input file
    //
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

    //
    // Parse input file
    //
    nlohmann::json j;

    input_file >> j;

    ExecutionStatistics stats;

    auto start = std::chrono::high_resolution_clock::now();
    std::unique_ptr<Mesh> mesh_ptr = j["mesh"];
    auto end = std::chrono::high_resolution_clock::now();
    stats.grid_building_time =
        std::chrono::duration<double>(end - start).count();

    vector<real_t> advection_x = j["advection_x"];
    vector<real_t> advection_y = j["advection_y"];
    vector<real_t> advection_z = j["advection_z"];
    vector<real_t> diffusion = j["diffusion"];
    vector<real_t> density = j["density"];

    real_t t_final = j["t_final"];
    real_t dt = j["dt"];

    start = std::chrono::high_resolution_clock::now();
    advection_diffusion_solver(
        *mesh_ptr, advection_x, advection_y, advection_z, diffusion, t_final,
        dt, density);
    end = std::chrono::high_resolution_clock::now();

    stats.trace_time = std::chrono::duration<double>(end - start).count();

    stats.total_time = stats.grid_building_time + stats.trace_time;

    //
    // Write output to file
    //
    nlohmann::json output;

    output["density"] = density;
    output["stats"] = stats;

    output_file << output.dump(4);

    return 0;
}