#ifndef NEPTUNE_IO_EXTRA_HPP
#define NEPTUNE_IO_EXTRA_HPP

#include "neptune/io/io_common.hpp"
#include "neptune/simulator/execution_statistics.hpp"
#include "neptune/simulator/simulationparams.hpp"

namespace neptune {

inline SimulationParams read_simulation_params(const H5::Group& params_group)
{
    SimulationParams params;
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
    if (H5Lexists(params_group.getId(), "algorithm", H5P_DEFAULT)) {
        std::string algorithm_str =
            read_string(params_group.openDataSet("algorithm"));
        if (algorithm_str == "MC") {
            params.algorithm = AlgorithmType::MC;
        }
        else if (algorithm_str == "KDMC") {
            params.algorithm = AlgorithmType::KDMC;
        }
        else if (algorithm_str == "KDMC_OLD") {
            params.algorithm = AlgorithmType::KDMC_OLD;
        }
        else if (algorithm_str == "ARCMC") {
            params.algorithm = AlgorithmType::ARCMC;
        }
        else {
            throw std::runtime_error(
                "Unknown algorithm type: " + algorithm_str);
        }
    }
    if(params.algorithm == AlgorithmType::KDMC || params.algorithm == AlgorithmType::KDMC_OLD) {
        if (H5Lexists(params_group.getId(), "kdmc_dt", H5P_DEFAULT)) {
            params.kdmc_dt = read_real_number(params_group.openDataSet("kdmc_dt"));
        }
        else {
            throw std::runtime_error(
                "Missing kdmc_dt parameter for KDMC algorithm");
        }
    }
    return params;
}

inline void write_execution_statistics(H5::H5File& output_file,
                                       const ExecutionStatistics& stats)
{
    H5::Group stats_group = output_file.createGroup("execution_statistics");
    H5::DataSet num_crossings_dataset = stats_group.createDataSet(
        "num_crossings", H5::PredType::NATIVE_INT64, H5::DataSpace());
    num_crossings_dataset.write(&stats.num_crossings,
                                H5::PredType::NATIVE_INT64);
    H5::DataSet num_collisions_dataset = stats_group.createDataSet(
        "num_collisions", H5::PredType::NATIVE_INT64, H5::DataSpace());
    num_collisions_dataset.write(&stats.num_collisions,
                                 H5::PredType::NATIVE_INT64);
    H5::DataSet total_particles_dataset = stats_group.createDataSet(
        "total_particles", H5::PredType::NATIVE_INT64, H5::DataSpace());
    total_particles_dataset.write(&stats.total_particles,
                                  H5::PredType::NATIVE_INT64);
    H5::DataSet grid_building_time_dataset = stats_group.createDataSet(
        "grid_building_time", H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
    grid_building_time_dataset.write(&stats.grid_building_time,
                                     H5::PredType::NATIVE_DOUBLE);
    H5::DataSet trace_time_dataset = stats_group.createDataSet(
        "trace_time", H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
    trace_time_dataset.write(&stats.trace_time, H5::PredType::NATIVE_DOUBLE);
    H5::DataSet total_time_dataset = stats_group.createDataSet(
        "total_time", H5::PredType::NATIVE_DOUBLE, H5::DataSpace());
    total_time_dataset.write(&stats.total_time, H5::PredType::NATIVE_DOUBLE);
    H5::DataSet num_kd_steps_dataset = stats_group.createDataSet(
        "num_kd_steps", H5::PredType::NATIVE_INT64, H5::DataSpace());
    num_kd_steps_dataset.write(&stats.num_kd_steps, H5::PredType::NATIVE_INT64);
    H5::DataSet num_kd_to_kinetic_switches_dataset =
        stats_group.createDataSet("num_kd_to_kinetic_switches",
                                  H5::PredType::NATIVE_INT64, H5::DataSpace());
    num_kd_to_kinetic_switches_dataset.write(&stats.num_kd_to_kinetic_switches,
                                             H5::PredType::NATIVE_INT64);
}
}  // namespace neptune

#endif