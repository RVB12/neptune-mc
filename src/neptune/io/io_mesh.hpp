#ifndef NEPTUNE_IO_MESH_HPP
#define NEPTUNE_IO_MESH_HPP

#include "neptune/geometry/mesh.hpp"
#include "neptune/geometry/nonuniform1Dmesh.hpp"
#include "neptune/geometry/nonuniformmesh.hpp"
#include "neptune/geometry/rotatednonuniformmesh.hpp"
#include "neptune/geometry/rotatedtrianglemesh.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/io/io_common.hpp"

namespace neptune {

/**
 * Read a vector of triangle indices from a 2D dataset. The dataset should be of
 * type int64_t and have 3 columns.
 */
inline vector<TriangleIndices> read_triangle_indices_vector(
    const H5::DataSet& dataset)
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
    else if (mesh_type == "non_uniform") {
        vector<real_t> cb_x = read_real_vector(mesh_group.openDataSet("cb_x"));
        vector<real_t> cb_y = read_real_vector(mesh_group.openDataSet("cb_y"));
        vector<real_t> cb_z = read_real_vector(mesh_group.openDataSet("cb_z"));
        NonUniformMesh mesh(cb_x, cb_y, cb_z);
        return std::make_unique<NonUniformMesh>(mesh);
    }
    else if (mesh_type == "non_uniform_1D") {
        vector<real_t> cb_x = read_real_vector(mesh_group.openDataSet("cb_x"));
        NonUniform1DMesh mesh(cb_x);
        return std::make_unique<NonUniform1DMesh>(mesh);
    }
    else if (mesh_type == "triangle") {
        vector<Vec2> vertices =
            read_vec2_vector(mesh_group.openDataSet("vertices"));
        vector<TriangleIndices> triangles =
            read_triangle_indices_vector(mesh_group.openDataSet("triangles"));
        if (mesh_group.exists("neighbors")) {
            vector<index_t> neighbors =
                read_index_vector(mesh_group.openDataSet("neighbors"));
            return std::make_unique<TriangleMesh>(
                TriangleMesh(vertices, triangles, neighbors));
        }
        else {
            return std::make_unique<TriangleMesh>(
                TriangleMesh(vertices, triangles));
        }
    }
    else if (mesh_type == "rotated_triangle") {
        vector<Vec2> vertices =
            read_vec2_vector(mesh_group.openDataSet("vertices"));
        vector<TriangleIndices> triangles =
            read_triangle_indices_vector(mesh_group.openDataSet("triangles"));
        if (mesh_group.exists("neighbors")) {
            vector<index_t> neighbors =
                read_index_vector(mesh_group.openDataSet("neighbors"));
            return std::make_unique<RotatedTriangleMesh>(
                RotatedTriangleMesh(vertices, triangles, neighbors));
        }
        else {
            return std::make_unique<RotatedTriangleMesh>(
                RotatedTriangleMesh(vertices, triangles));
        }
    }
    else if (mesh_type == "rotated_non_uniform") {
        vector<real_t> cb_r = read_real_vector(mesh_group.openDataSet("cb_r"));
        vector<real_t> cb_z = read_real_vector(mesh_group.openDataSet("cb_z"));
        RotatedNonUniformMesh mesh(cb_r, cb_z);
        return std::make_unique<RotatedNonUniformMesh>(mesh);
    }
    else {
        throw std::runtime_error("Unknown mesh type: " + mesh_type);
    }
}

}  // namespace neptune

#endif