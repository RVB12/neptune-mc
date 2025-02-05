#include <pybind11/pybind11.h>
#include <chrono>

#include <eigrene/simulator/simulator.hpp>
#include "eigrene/simulator/simulator_kdmc.hpp"
#include "eigrene/common/json_conversions.hpp"

#include <string>
#include <nlohmann/json.hpp>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

using namespace eigrene;

std::string simulate_json(std::string jsonstring){
    nlohmann::json j = nlohmann::json::parse(jsonstring);
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
    vector<real_t> density(mesh_ptr->n_cells());
    // TODO: read boundary conditions from json
    vector<BoundaryCondition> boundary_conditions = j["boundary_conditions"];



    start = std::chrono::high_resolution_clock::now();
    if(j.contains("dt")) {
        real_t dt = j["dt"];
        simulate_kdmc(*mesh_ptr, background, sources, density, boundary_conditions,
                      params, stats, dt);
    }
    else {
        simulate(*mesh_ptr, background, sources, density, boundary_conditions,
                      params, stats);
    }
    end = std::chrono::high_resolution_clock::now();

    stats.trace_time = std::chrono::duration<double>(end - start).count();
    stats.total_time = stats.grid_building_time + stats.trace_time;
    nlohmann::json output;
    output["density"] = density;
    output["stats"] = stats;
    return output.dump(4);
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        Eigrene - simple C++ monte carlo neutral code
        -----------------------

        .. currentmodule:: eigrene

        .. autosummary::
           :toctree: _generate
    )pbdoc";

    m.def("simulate", &simulate_json, R"pbdoc(
        Simulate particles and tally their paths
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
