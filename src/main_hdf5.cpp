#include <iostream>
#include <string>

#include "H5Cpp.h"
#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/collisions/simple_collisions.hpp"
#include "neptune/common/utils.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/background.hpp"

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
        throw std::runtime_error("Triangle mesh is not implemented yet");
    }
    else if (mesh_type == "rotated_triangle") {
        throw std::runtime_error(
            "Rotated triangle mesh is not implemented yet");
    }
    else if (mesh_type == "non_uniform") {
        throw std::runtime_error("Non-uniform mesh is not implemented yet");
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
                if (background_species == Species::Dplus &&
                    particle_species == Species::D) {
                    collisions.push_back(
                        std::make_unique<SimpleChargeExchangeCollision_D>());
                }
                else {
                    throw std::runtime_error(
                        "Unsupported species combination for simple charge "
                        "exchange collision: background species " +
                        speciestostring(background_species) +
                        ", particle species " +
                        speciestostring(particle_species));
                }
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
                if (background_species == Species::Dplus &&
                    particle_species == Species::D) {
                    collisions.push_back(
                        std::make_unique<SimpleIonisationCollision_D>());
                }
                else {
                    throw std::runtime_error(
                        "Unsupported species combination for simple ionisation "
                        "collision: background species " +
                        speciestostring(background_species) +
                        ", particle species " +
                        speciestostring(particle_species));
                }
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

    return 0;
}