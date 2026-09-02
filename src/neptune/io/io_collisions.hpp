#ifndef NEPTUNE_IO_COLLISIONS_HPP
#define NEPTUNE_IO_COLLISIONS_HPP

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/collisions/simple_collisions.hpp"
#include "neptune/collisions/user_defined_collisions.hpp"
#include "neptune/collisions/amjuel_collisions.hpp"
#include "neptune/io/io_common.hpp"

namespace neptune {

inline std::tuple<vector<std::unique_ptr<ChargeExchangeCollision>>,
                  vector<std::unique_ptr<IonisationCollision>>>
read_collisions(const H5::Group& collisions_group)
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
            else if (model == "user_defined") {
                Species background_species = read_species(
                    collision_group.openDataSet("background_species"));
                Species particle_species = read_species(
                    collision_group.openDataSet("particle_species"));
                vector<real_t> collision_rate = read_real_vector(
                    collision_group.openDataSet("collision_rate"));
                vector<real_t> velocity_mean_x = read_real_vector(
                    collision_group.openDataSet("velocity_mean_x"));
                vector<real_t> velocity_mean_y = read_real_vector(
                    collision_group.openDataSet("velocity_mean_y"));
                vector<real_t> velocity_mean_z = read_real_vector(
                    collision_group.openDataSet("velocity_mean_z"));
                vector<real_t> velocity_var = read_real_vector(
                    collision_group.openDataSet("velocity_var"));
                collisions.push_back(
                    std::make_unique<UserDefinedChargeExchangeCollision>(
                        background_species, particle_species, collision_rate,
                        velocity_mean_x, velocity_mean_y, velocity_mean_z,
                        velocity_var));
            }
            else if (model == "amjuel") {
                Species background_species = read_species(
                    collision_group.openDataSet("background_species"));
                Species particle_species = read_species(
                    collision_group.openDataSet("particle_species"));
                vector<real_t> coef_matrix = read_real_vector(
                    collision_group.openDataSet("coef_matrix"));
                collisions.push_back(
                    std::make_unique<AmjuelChargeExchangeCollision>(
                        background_species, particle_species, coef_matrix));
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
            else if (model == "user_defined") {
                Species background_species = read_species(
                    collision_group.openDataSet("background_species"));
                Species particle_species = read_species(
                    collision_group.openDataSet("particle_species"));
                vector<real_t> collision_rate = read_real_vector(
                    collision_group.openDataSet("collision_rate"));
                collisions.push_back(
                    std::make_unique<UserDefinedIonisationCollision>(
                        background_species, particle_species, collision_rate));
            }
            else if (model == "amjuel") {
                Species background_species = read_species(
                    collision_group.openDataSet("background_species"));
                Species particle_species = read_species(
                    collision_group.openDataSet("particle_species"));
                vector<real_t> coef_matrix = read_real_vector(
                    collision_group.openDataSet("coef_matrix"));
                collisions.push_back(
                    std::make_unique<AmjuelIonisationCollision>(
                        background_species, particle_species, coef_matrix));
            }
            else {
                throw std::runtime_error(
                    "Unsupported model for ionisation collision: " + model);
            }
        }

        if (H5Lexists(collision_group.getId(), "output_collision_info",
                      H5P_DEFAULT)) {
            bool output_collision_info =
                read_bool(collision_group.openDataSet("output_collision_info"));
            collisions[collisions.size() - 1]->output_collision_info =
                output_collision_info;
        }
    }

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

    return std::make_tuple(std::move(cx_collisions),
                           std::move(ionisation_collisions));
}

inline void write_collision_rates(
    const H5::H5File& output_file,
    const vector<std::unique_ptr<ChargeExchangeCollision>>& cx_collisions,
    const vector<std::unique_ptr<IonisationCollision>>& ionisation_collisions,
    const Mesh& mesh,
    const Background& background)
{
    bool output_collision_info = false;
    for (index_t i = 0; i < cx_collisions.size(); ++i) {
        if (cx_collisions[i]->output_collision_info) {
            output_collision_info = true;
            break;
        }
    }

    for (index_t i = 0; i < ionisation_collisions.size(); ++i) {
        if (ionisation_collisions[i]->output_collision_info) {
            output_collision_info = true;
            break;
        }
    }
    if (!output_collision_info) {
        return;
    }
    H5::Group collision_rates_group = output_file.createGroup("collision_info");
    for (index_t i = 0; i < cx_collisions.size(); ++i) {
        if (cx_collisions[i]->output_collision_info) {
            auto cx_group = collision_rates_group.createGroup(
                "charge_exchange_" + std::to_string(i));
            // Calculate the rates
            vector<real_t> rates(mesh.n_cells());
            vector<real_t> vel_mean_x(mesh.n_cells());
            vector<real_t> vel_mean_y(mesh.n_cells());
            vector<real_t> vel_mean_z(mesh.n_cells());
            vector<real_t> vel_var(mesh.n_cells());
            for (index_t j = 0; j < mesh.n_cells(); ++j) {
                // Make a particle with the background species and the cell
                // index
                Particle particle(0, 0, cx_collisions[i]->particle_species());
                particle.geometry.cell_id = j;
                Vec3 vel_mean = cx_collisions[i]->post_collision_velocity_mean(
                    particle, background);
                real_t vel_variance =
                    cx_collisions[i]->post_collision_velocity_var(particle,
                                                                  background);
                particle.geometry.velocity = vel_mean + Vec3(2*std::sqrt(vel_variance), 0, 0);
                real_t rate =
                    cx_collisions[i]->collision_rate(particle, background);
                vel_mean_x[j] = vel_mean.x;
                vel_mean_y[j] = vel_mean.y;
                vel_mean_z[j] = vel_mean.z;
                vel_var[j] = vel_variance;
                rates[j] = rate;
            }
            // Write the rates to the dataset
            hsize_t dim = rates.size();
            H5::DataSpace dataspace(1, &dim);  // 1D array
            H5::DataSet dataset = cx_group.createDataSet(
                "collision_rate", H5::PredType::NATIVE_DOUBLE, dataspace);
            dataset.write(rates.data(), H5::PredType::NATIVE_DOUBLE);
            H5::DataSet vel_mean_x_dataset =
                cx_group.createDataSet("post_collision_velocity_mean_x",
                                       H5::PredType::NATIVE_DOUBLE, dataspace);
            vel_mean_x_dataset.write(vel_mean_x.data(),
                                     H5::PredType::NATIVE_DOUBLE);
            H5::DataSet vel_mean_y_dataset =
                cx_group.createDataSet("post_collision_velocity_mean_y",
                                       H5::PredType::NATIVE_DOUBLE, dataspace);
            vel_mean_y_dataset.write(vel_mean_y.data(),
                                     H5::PredType::NATIVE_DOUBLE);
            H5::DataSet vel_mean_z_dataset =
                cx_group.createDataSet("post_collision_velocity_mean_z",
                                       H5::PredType::NATIVE_DOUBLE, dataspace);
            vel_mean_z_dataset.write(vel_mean_z.data(),
                                     H5::PredType::NATIVE_DOUBLE);
            // Write the velocity variance to the dataset
            H5::DataSet vel_var_dataset =
                cx_group.createDataSet("post_collision_velocity_variance",
                                       H5::PredType::NATIVE_DOUBLE, dataspace);
            vel_var_dataset.write(vel_var.data(), H5::PredType::NATIVE_DOUBLE);
        }
    }
    for (index_t i = 0; i < ionisation_collisions.size(); ++i) {
        if (ionisation_collisions[i]->output_collision_info) {
            auto ionisation_group = collision_rates_group.createGroup(
                "ionisation_" + std::to_string(i));
            // Calculate the rates
            vector<real_t> rates(mesh.n_cells());
            for (index_t j = 0; j < mesh.n_cells(); ++j) {
                // Make a particle with the background species and the cell
                // index
                Particle particle(0, 0,
                                  ionisation_collisions[i]->particle_species());
                particle.geometry.cell_id = j;
                Vec3 vel_mean = cx_collisions[i]->post_collision_velocity_mean(
                    particle, background);
                real_t vel_variance =
                    cx_collisions[i]->post_collision_velocity_var(particle,
                                                                  background);
                particle.geometry.velocity = vel_mean + Vec3(2*std::sqrt(vel_variance), 0, 0);
                real_t rate = ionisation_collisions[i]->collision_rate(
                    particle, background);
                rates[j] = rate;
            }
            // Write the rates to the dataset
            hsize_t dim = rates.size();
            H5::DataSpace dataspace(1, &dim);  // 1D array
            H5::DataSet dataset = ionisation_group.createDataSet(
                "collision_rate", H5::PredType::NATIVE_DOUBLE, dataspace);
            dataset.write(rates.data(), H5::PredType::NATIVE_DOUBLE);
        }
    }
}

}  // namespace neptune

#endif