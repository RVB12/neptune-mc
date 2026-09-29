#ifndef NEPTUNE_IO_COMMON_HPP
#define NEPTUNE_IO_COMMON_HPP

#include "H5Cpp.h"
#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/vec3.hpp"
#include "neptune/physics/species.hpp"

namespace neptune {

inline std::string get_dataset_path(const H5::DataSet& dataset)
{
    ssize_t size = dataset.getObjName(nullptr, 0);
    if (size < 0) return "<unknown>";

    std::string name(size + 1, '\0');
    dataset.getObjName(name.data(), size + 1);

    name.resize(size);  // remove trailing null
    return name;
}

/**
 * Read a real number from a scalar dataset.
 */
inline real_t read_real_number(const H5::DataSet& dataset)
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

/**
 * Read an integer number from a scalar dataset.
 */
inline index_t read_index_number(const H5::DataSet& dataset)
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

/**
 * Read a boolean value from a scalar dataset. The dataset should be of type
 * int64_t, where 0 is false and any non-zero value is true.
 */
inline bool read_bool(const H5::DataSet& dataset)
{
    index_t value = read_index_number(dataset);
    return value != 0;
}

/**
 * Read a vector of real numbers from a 1D dataset. The dataset should be of type double.
 */
inline vector<real_t> read_real_vector(const H5::DataSet& dataset)
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

/**
 * Read a vector of Vec2 objects from a 2D dataset. The dataset should be of type double and have 2 columns.
 */
inline vector<Vec2> read_vec2_vector(const H5::DataSet& dataset)
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

/**
 * Read a vector of index_t from a 1D dataset. The dataset should be of type int64_t.
 */
inline vector<index_t> read_index_vector(const H5::DataSet& dataset)
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

/**
 * Read a string from a scalar dataset. The dataset should be of type string.
 */
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

/**
 * Read a Species from a scalar dataset. The dataset should be of type string.
 */
inline Species read_species(const H5::DataSet& dataset)
{
    std::string species_str = read_string(dataset);
    return stringtospecies(species_str);
}

}  // namespace neptune

#endif