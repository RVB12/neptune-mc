#ifndef NEPTUNE_IO_BACKGROUND_HPP
#define NEPTUNE_IO_BACKGROUND_HPP

#include "neptune/io/io_common.hpp"
#include "neptune/physics/background.hpp"

namespace neptune {

inline Background read_background(const H5::Group& background_group)
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

}  // namespace neptune

#endif