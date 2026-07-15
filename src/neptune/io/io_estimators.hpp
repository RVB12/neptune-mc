#ifndef NEPTUNE_IO_ESTIMATORS_HPP
#define NEPTUNE_IO_ESTIMATORS_HPP

#include "neptune/estimation/estimator.hpp"
#include "neptune/estimation/estimatorgroup.hpp"
#include "neptune/io/io_common.hpp"

namespace neptune {

inline EstimatorGroup read_estimators(const H5::Group& estimators_group,
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

inline void write_estimators(H5::H5File& output_file,
                             const EstimatorGroup& estimator_group)
{
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
        bool has_kdmc =
            !estimator_group
                 .kdmc_estimators_by_species[static_cast<index_t>(species)]
                 .empty();
        if (!has_tl && !has_final && !has_kdmc) {
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

        if (has_kdmc) {
            H5::Group kdmc_group = species_group.createGroup("kdmc_fluid");

            for (auto estimator : estimator_group.kdmc_estimators_by_species
                                      [static_cast<index_t>(species)]) {
                std::string dataset_name = QOItostring(estimator->qoi);

                vector<real_t> estimate_vector = estimator->get_estimate();

                hsize_t dim = estimate_vector.size();
                H5::DataSpace dataspace(1, &dim);  // 1D array

                H5::DataSet dataset = kdmc_group.createDataSet(
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
                    H5::DataSet variance_dataset = kdmc_group.createDataSet(
                        dataset_name + "_variance", H5::PredType::NATIVE_DOUBLE,
                        dataspace_var);

                    vector<real_t> variance_vector = estimator->get_variance();

                    variance_dataset.write(variance_vector.data(),
                                           H5::PredType::NATIVE_DOUBLE);
                }
            }
        }
    }
}

}  // namespace neptune

#endif