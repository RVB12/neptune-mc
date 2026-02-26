#include <chrono>
#include <fstream>

#include "neptune/common/version.hpp"
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
    if (argc == 2) {
        std::string argument = argv[1];
        if (argument == "-v" || argument == "--version") {
            std::cout << "Neptune version: " << NEPTUNE_VERSION
                      << std::endl;
            return 0;
        }
    }
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
    vector<real_t> diffusion_rate = j["diffusion_rate"];
    vector<real_t> diffusion_factor = j["diffusion_factor"];
    vector<real_t> density = j["density"];

    real_t t_final = j["t_final"];
    real_t dt = j["dt"];

    vector<real_t> density_integral(density.size(), 0.0);
    vector<real_t> momentum_x_integral(density.size(), 0.0);
    vector<real_t> momentum_y_integral(density.size(), 0.0);
    vector<real_t> momentum_z_integral(density.size(), 0.0);
    vector<real_t> energy_integral(density.size(), 0.0);
    vector<real_t> energy_x_integral(density.size(), 0.0);

    vector<std::unique_ptr<FluidBoundaryCondition>> boundary_conditions =
        j["boundary_conditions"];

    start = std::chrono::high_resolution_clock::now();
    advection_diffusion_solver(*mesh_ptr, boundary_conditions, advection_x,
                               advection_y, advection_z, diffusion_rate,
                               diffusion_factor, t_final, dt, density,
                               density_integral, momentum_x_integral, momentum_y_integral,
                               momentum_z_integral, energy_integral, energy_x_integral);
    end = std::chrono::high_resolution_clock::now();

    stats.trace_time = std::chrono::duration<double>(end - start).count();

    stats.total_time = stats.grid_building_time + stats.trace_time;

    //
    // Write output to file
    //
    nlohmann::json output;

    output["density"] = density;
    output["density_integral"] = density_integral;
    output["momentum_x_integral"] = momentum_x_integral;
    output["momentum_y_integral"] = momentum_y_integral;
    output["momentum_z_integral"] = momentum_z_integral;
    output["energy_integral"] = energy_integral;
    output["energy_x_integral"] = energy_x_integral;
    output["stats"] = stats;

    output_file << output.dump(4);

    return 0;
}