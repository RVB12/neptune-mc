#include <iostream>
#include <string>

#include "H5Cpp.h"
#include "neptune/boundaryconditions/boundarycondition.hpp"
#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/collisions/simple_collisions.hpp"
#include "neptune/common/utils.hpp"
#include "neptune/estimation/estimator.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/rotatedtrianglemesh.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/nonuniformmesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/background.hpp"
#include "neptune/physics/species.hpp"
#include "neptune/simulator/simulator.hpp"
#include "neptune/sources/recombinationsource.hpp"
#include "neptune/sources/source.hpp"
#include "neptune/sources/volumesource.hpp"
#include "neptune/sources/wallreflectionsource.hpp"
#include "neptune/boundaryconditions/trim.hpp"

using namespace neptune;

std::string get_dataset_path(const H5::DataSet& dataset)
{
    ssize_t size = dataset.getObjName(nullptr, 0);
    if (size < 0) return "<unknown>";

    std::string name(size + 1, '\0');
    dataset.getObjName(name.data(), size + 1);

    name.resize(size);  // remove trailing null
    return name;
}

real_t read_real_number(const H5::DataSet& dataset)
{
    H5::DataSpace dataspace = dataset.getSpace();
    if (dataspace.getSimpleExtentNdims() != 0) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not scalar");
    }
    if (dataset.getDataType() != H5::PredType::NATIVE_DOUBLE) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not of type double");
    }
    real_t value;
    dataset.read(&value, H5::PredType::NATIVE_DOUBLE);
    return value;
}

index_t read_index_number(const H5::DataSet& dataset)
{
    H5::DataSpace dataspace = dataset.getSpace();
    if (dataspace.getSimpleExtentNdims() != 0) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not scalar");
    }
    if (dataset.getDataType() != H5::PredType::NATIVE_INT64) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not of type int64_t");
    }
    index_t value;
    dataset.read(&value, H5::PredType::NATIVE_INT64);
    return value;
}

bool read_bool(const H5::DataSet& dataset)
{
    index_t value = read_index_number(dataset);
    return value != 0;
}

vector<real_t> read_real_vector(const H5::DataSet& dataset)
{
    H5::DataSpace dataspace = dataset.getSpace();
    if (dataspace.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not 1D");
    }
    // check if the dataset is of type double
    if (dataset.getDataType() != H5::PredType::NATIVE_DOUBLE) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not of type double");
    }
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, nullptr);
    vector<real_t> data(dims[0]);
    dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    return data;
}

vector<Vec2> read_vec2_vector(const H5::DataSet& dataset)
{
    H5::DataSpace dataspace = dataset.getSpace();
    if (dataspace.getSimpleExtentNdims() != 2) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not 2D");
    }
    // check if the dataset is of type double
    if (dataset.getDataType() != H5::PredType::NATIVE_DOUBLE) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not of type double");
    }
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims, nullptr);
    if (dims[1] != 2) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " does not have 2 columns");
    }
    vector<Vec2> data(dims[0]);
    vector<real_t> raw_data(dims[0] * dims[1]);
    dataset.read(raw_data.data(), H5::PredType::NATIVE_DOUBLE);
    for (index_t i = 0; i < dims[0]; ++i) {
        data[i] = Vec2(raw_data[2 * i], raw_data[2 * i + 1]);
    }
    return data;
}

vector<TriangleIndices> read_triangle_indices_vector(const H5::DataSet& dataset)
{
    H5::DataSpace dataspace = dataset.getSpace();
    if (dataspace.getSimpleExtentNdims() != 2) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not 2D");
    }
    // check if the dataset is of type double
    if (dataset.getDataType() != H5::PredType::NATIVE_INT64) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not of type int64_t");
    }
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims, nullptr);
    if (dims[1] != 3) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " does not have 3 columns");
    }
    vector<TriangleIndices> data(dims[0]);
    vector<index_t> raw_data(dims[0] * dims[1]);
    dataset.read(raw_data.data(), H5::PredType::NATIVE_INT64);
    for (index_t i = 0; i < dims[0]; ++i) {
        data[i] = TriangleIndices(raw_data[3 * i], raw_data[3 * i + 1],
                                  raw_data[3 * i + 2]);
    }
    return data;
}

vector<index_t> read_index_vector(const H5::DataSet& dataset)
{
    H5::DataSpace dataspace = dataset.getSpace();
    if (dataspace.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not 1D");
    }
    // check if the dataset is of type int64_t
    if (dataset.getDataType() != H5::PredType::NATIVE_INT64) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not of type int64_t");
    }
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, nullptr);
    vector<index_t> data(dims[0]);
    dataset.read(data.data(), H5::PredType::NATIVE_INT64);
    return data;
}

std::string read_string(const H5::DataSet& dataset)
{
    H5::DataSpace dataspace = dataset.getSpace();
    if (dataspace.getSimpleExtentNdims() != 0) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not scalar");
    }

    H5::DataType dtype = dataset.getDataType();
    if (dtype.getClass() != H5T_STRING) {
        throw std::runtime_error("Dataset: " + get_dataset_path(dataset) +
                                 " is not a string");
    }

    H5::StrType strType = dataset.getStrType();

    if (strType.isVariableStr()) {
        char* buffer = nullptr;
        dataset.read(&buffer, strType);

        std::string result(buffer);
        free(buffer);

        return result;
    }
    else {
        size_t size = strType.getSize();

        std::string buffer(size, '\0');
        dataset.read(buffer.data(), strType);

        // Remove trailing nulls
        buffer.erase(buffer.find('\0'));

        return buffer;
    }
}

Species read_species(const H5::DataSet& dataset)
{
    std::string species_str = read_string(dataset);
    return stringtospecies(species_str);
}

Background read_background(const H5::Group& background_group)
{
    Background background;
    for (index_t species_index = 0; species_index < NUM_SPECIES;
         species_index++) {
        Species species = static_cast<Species>(species_index);
        std::string species_name = speciestostring(species);
        // Check if the species group exists
        if (H5Lexists(background_group.getId(), species_name.c_str(),
                      H5P_DEFAULT)) {
            H5::Group species_group = background_group.openGroup(species_name);

            H5::DataSet density_dataset = species_group.openDataSet("density");
            H5::DataSet velocity_x_dataset =
                species_group.openDataSet("velocity_x");
            H5::DataSet velocity_y_dataset =
                species_group.openDataSet("velocity_y");
            H5::DataSet velocity_z_dataset =
                species_group.openDataSet("velocity_z");
            H5::DataSet temperature_dataset =
                species_group.openDataSet("temperature");
            vector<real_t> density = read_real_vector(density_dataset);
            vector<real_t> velocity_x = read_real_vector(velocity_x_dataset);
            vector<real_t> velocity_y = read_real_vector(velocity_y_dataset);
            vector<real_t> velocity_z = read_real_vector(velocity_z_dataset);
            vector<real_t> temperature = read_real_vector(temperature_dataset);

            background.set_species_data(
                species, std::move(density), std::move(velocity_x),
                std::move(velocity_y), std::move(velocity_z),
                std::move(temperature));
        }
    }
    return background;
}

std::unique_ptr<Mesh> read_mesh(const H5::Group& mesh_group)
{
    // First, read the mesh type from the attribute "type" of the mesh group
    std::string mesh_type = read_string(mesh_group.openDataSet("type"));
    if (mesh_type == "uniform") {
        real_t Lx = read_real_number(mesh_group.openDataSet("Lx"));
        real_t Ly = read_real_number(mesh_group.openDataSet("Ly"));
        real_t Lz = read_real_number(mesh_group.openDataSet("Lz"));
        index_t nx = read_index_number(mesh_group.openDataSet("nx"));
        index_t ny = read_index_number(mesh_group.openDataSet("ny"));
        index_t nz = read_index_number(mesh_group.openDataSet("nz"));
        VoxelMesh mesh(Lx, Ly, Lz, nx, ny, nz);
        return std::make_unique<VoxelMesh>(mesh);
    }
    else if (mesh_type == "triangle") {
        vector<Vec2> vertices =
            read_vec2_vector(mesh_group.openDataSet("vertices"));
        vector<TriangleIndices> triangles =
            read_triangle_indices_vector(mesh_group.openDataSet("triangles"));
        TriangleMesh mesh(vertices, triangles);
        return std::make_unique<TriangleMesh>(mesh);
    }
    else if (mesh_type == "rotated_triangle") {
        vector<Vec2> vertices =
            read_vec2_vector(mesh_group.openDataSet("vertices"));
        vector<TriangleIndices> triangles =
            read_triangle_indices_vector(mesh_group.openDataSet("triangles"));
        RotatedTriangleMesh mesh(vertices, triangles);
        return std::make_unique<RotatedTriangleMesh>(mesh);
    }
    else if (mesh_type == "non_uniform") {
        vector<real_t> cb_x =
            read_real_vector(mesh_group.openDataSet("cb_x"));
        vector<real_t> cb_y =
            read_real_vector(mesh_group.openDataSet("cb_y"));
        vector<real_t> cb_z =
            read_real_vector(mesh_group.openDataSet("cb_z"));
        NonUniformMesh mesh(cb_x, cb_y, cb_z);
        return std::make_unique<NonUniformMesh>(mesh);
    }
    else if (mesh_type == "rotated_non_uniform") {
        throw std::runtime_error(
            "Rotated non-uniform mesh is not implemented yet");
    }
    else {
        throw std::runtime_error("Unknown mesh type: " + mesh_type);
    }
}

vector<std::unique_ptr<Collision>> read_collisions(
    const H5::Group& collisions_group)
{
    std::vector<std::unique_ptr<Collision>> collisions;

    // Each subgroup of the collisions group corresponds to a collision process
    for (hsize_t i = 0; i < collisions_group.getNumObjs(); i++) {
        std::string collision_name = collisions_group.getObjnameByIdx(i);
        H5::Group collision_group = collisions_group.openGroup(collision_name);
        std::string collision_type =
            read_string(collision_group.openDataSet("type"));

        if (collision_type == "charge_exchange") {
            // Read charge exchange collision data
            std::string model =
                read_string(collision_group.openDataSet("model"));
            if (model == "simple") {
                Species background_species = read_species(
                    collision_group.openDataSet("background_species"));
                Species particle_species = read_species(
                    collision_group.openDataSet("particle_species"));
                collisions.push_back(
                    std::make_unique<SimpleChargeExchangeCollision>(
                        background_species, particle_species));
            }
            else {
                throw std::runtime_error(
                    "Unsupported model for charge exchange collision: " +
                    model);
            }
        }
        else if (collision_type == "ionisation") {
            // Read ionisation collision data
            std::string model =
                read_string(collision_group.openDataSet("model"));
            if (model == "simple") {
                Species background_species = read_species(
                    collision_group.openDataSet("background_species"));
                Species particle_species = read_species(
                    collision_group.openDataSet("particle_species"));
                collisions.push_back(
                    std::make_unique<SimpleIonisationCollision>(
                        background_species, particle_species));
            }
            else {
                throw std::runtime_error(
                    "Unsupported model for ionisation collision: " + model);
            }
        }
    }

    // Convert std::vector to neptune::vector
    vector<std::unique_ptr<Collision>> collisions_vector(collisions.size());
    for (size_t i = 0; i < collisions.size(); i++) {
        collisions_vector[i] = std::move(collisions[i]);
    }

    return collisions_vector;
}

vector<std::unique_ptr<Source>> read_sources(const H5::Group& sources_group,
                                             const Mesh& mesh,
                                             const Background& background)
{
    std::vector<std::unique_ptr<Source>> sources;

    // Go over all the subgroups of the sources group, each subgroup corresponds
    // to a source
    for (hsize_t i = 0; i < sources_group.getNumObjs(); i++) {
        if (!sources_group.getObjTypeByIdx(i) == H5G_GROUP) {
            continue;
        }
        std::string source_name = sources_group.getObjnameByIdx(i);
        H5::Group source_group = sources_group.openGroup(source_name);
        std::string source_type = read_string(source_group.openDataSet("type"));

        if (source_type == "recombination") {
            // Read recombination source data
            Species background_species =
                read_species(source_group.openDataSet("background_species"));
            Species particle_species =
                read_species(source_group.openDataSet("particle_species"));
            vector<index_t> cell_indices =
                read_index_vector(source_group.openDataSet("cell_indices"));
            for (index_t cell_index : cell_indices) {
                sources.push_back(std::make_unique<RecombinationSource>(
                    0, cell_index, background_species, particle_species, mesh,
                    background));
            }
        }
        else if (source_type == "wall_recycling") {
            // Read wall recycling source data
            vector<index_t> cell_indices =
                read_index_vector(source_group.openDataSet("cell_indices"));
            vector<index_t> surface_indices =
                read_index_vector(source_group.openDataSet("surface_indices"));
            real_t sheath_potential =
                read_real_number(source_group.openDataSet("sheath_potential"));
            real_t franck_condon_temp = read_real_number(
                source_group.openDataSet("franck_condon_temp"));
            real_t Rt = read_real_number(source_group.openDataSet("Rt"));
            real_t Rr = read_real_number(source_group.openDataSet("Rr"));
            Species background_species =
                read_species(source_group.openDataSet("background_species"));
            Species particle_species =
                read_species(source_group.openDataSet("particle_species"));
            for (size_t j = 0; j < cell_indices.size(); j++) {
                index_t cell_index = cell_indices[j];
                index_t surface_index = surface_indices[j];
                sources.push_back(std::make_unique<WallReflectionSource>(
                    0, cell_index, surface_index, sheath_potential,
                    franck_condon_temp, Rt, Rr, background_species,
                    particle_species, mesh, background));
            }
        }
        else if (source_type == "volume") {
            Species species = read_species(source_group.openDataSet("species"));
            vector<index_t> cell_indices =
                read_index_vector(source_group.openDataSet("cell_indices"));
            vector<real_t> intensities =
                read_real_vector(source_group.openDataSet("intensities"));
            vector<real_t> vel_x =
                read_real_vector(source_group.openDataSet("velocity_x"));
            vector<real_t> vel_y =
                read_real_vector(source_group.openDataSet("velocity_y"));
            vector<real_t> vel_z =
                read_real_vector(source_group.openDataSet("velocity_z"));
            vector<real_t> vel_std =
                read_real_vector(source_group.openDataSet("velocity_std"));
            for (size_t j = 0; j < cell_indices.size(); j++) {
                index_t cell_index = cell_indices[j];
                real_t intensity = intensities[j];
                real_t velocity_x = vel_x[j];
                real_t velocity_y = vel_y[j];
                real_t velocity_z = vel_z[j];
                real_t velocity_std = vel_std[j];
                sources.push_back(std::make_unique<VolumeSource>(
                    0, species, intensity, cell_index,
                    Vec3(velocity_x, velocity_y, velocity_z), velocity_std));
            }
        }
        else {
            throw std::runtime_error("Unknown source type: " + source_type);
        }
    }

    // Automatic stratification, assign number of particles according to the
    // intensity of the sources
    index_t n_particles =
        read_index_number(sources_group.openDataSet("n_particles"));
    // real_t total_intensity = 0;
    // for (const auto& source : sources) {
    //     total_intensity += source->intensity();
    // }
    // for (auto& source : sources) {
    //     source->num_particles_ = static_cast<index_t>(
    //         n_particles * source->intensity() / total_intensity);
    // }
    for (auto& source : sources) {
        source->num_particles_ = static_cast<index_t>(
            n_particles / sources.size());  // equal number of particles for each source
    }

    // Convert std::vector to neptune::vector
    vector<std::unique_ptr<Source>> sources_vector(sources.size());
    for (size_t i = 0; i < sources.size(); i++) {
        sources_vector[i] = std::move(sources[i]);
    }
    return sources_vector;
}

EstimatorGroup read_estimators(const H5::Group& estimators_group,
                               index_t n_cells)
{
    std::vector<std::tuple<Species, QOI, EstimatorType, bool>>
        estimators_definitions;

    // Go over all the subgroups of the estimators group, each subgroup
    // corresponds to an estimator
    for (hsize_t i = 0; i < estimators_group.getNumObjs(); i++) {
        if (!estimators_group.getObjTypeByIdx(i) == H5G_GROUP) {
            continue;
        }
        std::string estimator_name = estimators_group.getObjnameByIdx(i);
        H5::Group estimator_group = estimators_group.openGroup(estimator_name);
        Species species = read_species(estimator_group.openDataSet("species"));
        QOI qoi = stringtoQOI(read_string(estimator_group.openDataSet("qoi")));
        EstimatorType estimator_type = stringtoEstimatorType(
            read_string(estimator_group.openDataSet("estimator_type")));
        bool estimate_variance = false;
        if (H5Lexists(estimator_group.getId(), "estimate_variance",
                      H5P_DEFAULT)) {
            bool estimate_variance =
                read_bool(estimator_group.openDataSet("estimate_variance"));
        }

        estimators_definitions.emplace_back(species, qoi, estimator_type,
                                            estimate_variance);
    }

    return EstimatorGroup(estimators_definitions, n_cells);
}

vector<std::unique_ptr<BoundaryCondition>> read_boundary_conditions(
    const H5::Group& boundary_conditions_group)
{
    std::vector<std::unique_ptr<BoundaryCondition>> boundary_conditions;

    // Go over all the subgroups of the boundary conditions group, each subgroup
    // corresponds to a boundary condition
    for (hsize_t i = 0; i < boundary_conditions_group.getNumObjs(); i++) {
        if (!boundary_conditions_group.getObjTypeByIdx(i) == H5G_GROUP) {
            continue;
        }
        std::string bc_name = boundary_conditions_group.getObjnameByIdx(i);
        H5::Group bc_group = boundary_conditions_group.openGroup(bc_name);
        std::string bc_type = read_string(bc_group.openDataSet("type"));

        if (bc_type == "specular") {
            boundary_conditions.push_back(std::make_unique<SpecularBC>());
        }
        else if (bc_type == "absorbing") {
            boundary_conditions.push_back(std::make_unique<AbsorbingBC>());
        }
        else if (bc_type == "partially_specular") {
            real_t absorption_coefficient = read_real_number(
                bc_group.openDataSet("absorption_coefficient"));
            boundary_conditions.push_back(
                std::make_unique<PartiallySpecularBC>(absorption_coefficient));
        }
        else {
            throw std::runtime_error("Unknown boundary condition type: " +
                                     bc_type);
        }
    }

    // Convert std::vector to neptune::vector
    vector<std::unique_ptr<BoundaryCondition>> bc_vector(
        boundary_conditions.size());
    for (size_t i = 0; i < boundary_conditions.size(); i++) {
        bc_vector[i] = std::move(boundary_conditions[i]);
    }
    return bc_vector;
}

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

    vector<std::unique_ptr<Collision>> collisions =
        read_collisions(input_file.openGroup("/collisions"));

    // Split collisions into CX and ionisation
    index_t n_cx_collisions = 0;
    index_t n_ionisation_collisions = 0;
    for (const auto& collision : collisions) {
        if (dynamic_cast<ChargeExchangeCollision*>(collision.get()) !=
            nullptr) {
            n_cx_collisions++;
        }
        else if (dynamic_cast<IonisationCollision*>(collision.get()) !=
                 nullptr) {
            n_ionisation_collisions++;
        }
    }
    vector<std::unique_ptr<ChargeExchangeCollision>> cx_collisions(
        n_cx_collisions);
    vector<std::unique_ptr<IonisationCollision>> ionisation_collisions(
        n_ionisation_collisions);
    index_t cx_index = 0;
    index_t ionisation_index = 0;
    for (auto& collision : collisions) {
        if (auto* cx_collision =
                dynamic_cast<ChargeExchangeCollision*>(collision.get())) {
            cx_collisions[cx_index++] =
                std::unique_ptr<ChargeExchangeCollision>(cx_collision);
            collision.release();  // release ownership from the original vector
        }
        else if (auto* ionisation_collision =
                     dynamic_cast<IonisationCollision*>(collision.get())) {
            ionisation_collisions[ionisation_index++] =
                std::unique_ptr<IonisationCollision>(ionisation_collision);
            collision.release();  // release ownership from the original vector
        }
    }

    vector<std::unique_ptr<Source>> sources =
        read_sources(input_file.openGroup("/sources"), *mesh, background);

    EstimatorGroup estimator_group =
        read_estimators(input_file.openGroup("/estimators"), mesh->n_cells());

    // Read TRIM database from the databases/TRIM folder in the project directory
    // This is a hardcoded path until I can figure out a better solution
    TrimReflectionData trim_data = read_trim_data("/home/thijs/Documents/Git/Monte-Carlo/neptune/databases/TRIM/D_on_W");

    vector<std::unique_ptr<BoundaryCondition>> boundary_conditions =
        read_boundary_conditions(input_file.openGroup("/boundary_conditions"));

    SimulationParams params;
    H5::Group params_group = input_file.openGroup("/params");
    if (H5Lexists(params_group.getId(), "max_steps", H5P_DEFAULT)) {
        params.max_steps =
            read_index_number(params_group.openDataSet("max_steps"));
    }
    if (H5Lexists(params_group.getId(), "seed", H5P_DEFAULT)) {
        params.seed = read_index_number(params_group.openDataSet("seed"));
    }
    if (H5Lexists(params_group.getId(), "t_final", H5P_DEFAULT)) {
        params.t_final = read_real_number(params_group.openDataSet("t_final"));
    }
    if (H5Lexists(params_group.getId(), "n_batches", H5P_DEFAULT)) {
        params.n_batches =
            read_index_number(params_group.openDataSet("n_batches"));
    }

    ExecutionStatistics stats;

    simulate(*mesh, background, cx_collisions, ionisation_collisions, sources,
             estimator_group, boundary_conditions, params, stats);

    // Now write the results to the output file
    H5::H5File output_file(output_file_name, H5F_ACC_TRUNC);

    for (const auto species : {Species::H, Species::D}) {
        bool has_tl =
            !estimator_group
                 .track_length_estimators_by_species[static_cast<index_t>(
                     species)]
                 .empty();
        bool has_final =
            !estimator_group
                 .final_estimators_by_species[static_cast<index_t>(species)]
                 .empty();
        if (!has_tl && !has_final) {
            continue;  // skip species for which we have no estimators
        }

        std::string species_name = speciestostring(species);
        H5::Group species_group = output_file.createGroup(species_name);

        if (has_tl) {
            H5::Group tl_group = species_group.createGroup("track_length");

            for (auto estimator :
                 estimator_group
                     .track_length_estimators_by_species[static_cast<index_t>(
                         species)]) {
                std::string dataset_name = QOItostring(estimator->qoi);

                vector<real_t> estimate_vector = estimator->get_estimate();

                hsize_t dim = estimate_vector.size();
                H5::DataSpace dataspace(1, &dim);  // 1D array

                H5::DataSet dataset = tl_group.createDataSet(
                    dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);

                dataset.write(estimate_vector.data(),
                              H5::PredType::NATIVE_DOUBLE);

                // Add unit attribute to the dataset
                std::string unit = getUnit(estimator->qoi);
                H5::StrType strType(H5::PredType::C_S1, unit.size());
                H5::Attribute unit_attr =
                    dataset.createAttribute("unit", strType, H5::DataSpace());
                unit_attr.write(strType, unit);

                if (estimator->estimate_variance) {
                    H5::DataSpace dataspace_var(1, &dim);  // 1D array
                    H5::DataSet variance_dataset = tl_group.createDataSet(
                        dataset_name + "_variance", H5::PredType::NATIVE_DOUBLE,
                        dataspace_var);

                    vector<real_t> variance_vector = estimator->get_variance();

                    variance_dataset.write(variance_vector.data(),
                                           H5::PredType::NATIVE_DOUBLE);
                }
            }
        }

        if (has_final) {
            H5::Group final_group = species_group.createGroup("final_time");

            for (auto estimator : estimator_group.final_estimators_by_species
                                      [static_cast<index_t>(species)]) {
                std::string dataset_name = QOItostring(estimator->qoi);

                vector<real_t> estimate_vector = estimator->get_estimate();

                hsize_t dim = estimate_vector.size();
                H5::DataSpace dataspace(1, &dim);  // 1D array

                H5::DataSet dataset = final_group.createDataSet(
                    dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);

                dataset.write(estimate_vector.data(),
                              H5::PredType::NATIVE_DOUBLE);

                // Add unit attribute to the dataset
                std::string unit = getUnit(estimator->qoi);
                H5::StrType strType(H5::PredType::C_S1, unit.size());
                H5::Attribute unit_attr =
                    dataset.createAttribute("unit", strType, H5::DataSpace());
                unit_attr.write(strType, unit);

                if (estimator->estimate_variance) {
                    H5::DataSpace dataspace_var(1, &dim);  // 1D array
                    H5::DataSet variance_dataset = final_group.createDataSet(
                        dataset_name + "_variance", H5::PredType::NATIVE_DOUBLE,
                        dataspace_var);

                    vector<real_t> variance_vector = estimator->get_variance();

                    variance_dataset.write(variance_vector.data(),
                                           H5::PredType::NATIVE_DOUBLE);
                }
            }
        }
    }

    return 0;
}