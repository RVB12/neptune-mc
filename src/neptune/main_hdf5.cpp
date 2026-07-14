#include <iostream>
#include <string>

#include "neptune/boundaryconditions/boundarycondition.hpp"
#include "neptune/boundaryconditions/trim.hpp"
#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/collisions/simple_collisions.hpp"
#include "neptune/collisions/user_defined_collisions.hpp"
#include "neptune/common/utils.hpp"
#include "neptune/estimation/estimator.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/nonuniformmesh.hpp"
#include "neptune/geometry/rotatedtrianglemesh.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/io/io_background.hpp"
#include "neptune/io/io_boundaryconditions.hpp"
#include "neptune/io/io_collisions.hpp"
#include "neptune/io/io_common.hpp"
#include "neptune/io/io_estimators.hpp"
#include "neptune/io/io_mesh.hpp"
#include "neptune/io/io_sources.hpp"
#include "neptune/io/io_extra.hpp"
#include "neptune/physics/background.hpp"
#include "neptune/physics/species.hpp"
#include "neptune/simulator/simulator.hpp"
#include "neptune/sources/recombinationsource.hpp"
#include "neptune/sources/source.hpp"
#include "neptune/sources/volumesource.hpp"
#include "neptune/sources/wallreflectionsource.hpp"

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
            std::cout << "Neptune version: " << NEPTUNE_VERSION << std::endl;
            return 0;
        }
    }
    if (argc != 3) {
        std::cout << "No input file specified, using default input.h5 and "
                     "output.h5"
                  << std::endl;
        input_file_name = "input.h5";
        output_file_name = "output.h5";
    }
    else {
        input_file_name = argv[1];
        output_file_name = argv[2];
    }

    H5::H5File input_file(input_file_name, H5F_ACC_RDONLY);

    Background background =
        read_background(input_file.openGroup("/background"));

    std::unique_ptr<Mesh> mesh = read_mesh(input_file.openGroup("/mesh"));

    auto [cx_collisions, ionisation_collisions] =
        read_collisions(input_file.openGroup("/collisions"));

    vector<std::unique_ptr<Source>> sources =
        read_sources(input_file.openGroup("/sources"), *mesh, background);

    EstimatorGroup estimator_group =
        read_estimators(input_file.openGroup("/estimators"), mesh->n_cells());

    // Read TRIM database from the databases/TRIM folder in the project
    // directory This is a hardcoded path until I can figure out a better
    // solution TrimReflectionData trim_data =
    // read_trim_data("/home/thijs/Documents/Git/Monte-Carlo/neptune/databases/TRIM/D_on_W");

    vector<std::unique_ptr<BoundaryCondition>> boundary_conditions =
        read_boundary_conditions(input_file.openGroup("/boundary_conditions"));

    SimulationParams params =
        read_simulation_params(input_file.openGroup("/params"));

    ExecutionStatistics stats;
    simulate(*mesh, background, cx_collisions, ionisation_collisions, sources,
             estimator_group, boundary_conditions, params, stats);

    // Now write the results to the output file
    H5::H5File output_file(output_file_name, H5F_ACC_TRUNC);

    //
    // Output execution statistics
    //
    write_execution_statistics(output_file, stats);

    // Output estimators
    write_estimators(output_file, estimator_group);

    // Output collision rates
    // (this is optional for debugging/post-processing)
    write_collision_rates(output_file, cx_collisions, ionisation_collisions,
                          *mesh, background);

    return 0;
}