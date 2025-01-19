#include <pybind11/pybind11.h>

#include <eiron2/simulator/simulator.hpp>
#include "eiron2/common/json_conversions.hpp"

#include <string>
#include <nlohmann/json.hpp>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

using namespace eiron2;

std::string simulate_json(std::string jsonstring){
    nlohmann::json j = nlohmann::json::parse(jsonstring);
    std::unique_ptr<Mesh> mesh_ptr = j["mesh"];
    PlasmaBackground background = j["background"];
    vector<std::unique_ptr<Source>> sources = j["sources"];
    SimulationParams params;
    params.max_steps = j["max_steps"];
    params.seed = j["seed"];
    vector<real_t> density(mesh_ptr->n_cells());
    // TODO: read boundary conditions from json
    vector<BoundaryCondition> boundary_conditions = j["boundary_conditions"];
    simulate(*mesh_ptr, background, sources, density, boundary_conditions, params);
    nlohmann::json output;
    output["density"] = density;
    return output.dump(4);
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        Eiron2 - simple C++ monte carlo neutral code
        -----------------------

        .. currentmodule:: eiron2

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
