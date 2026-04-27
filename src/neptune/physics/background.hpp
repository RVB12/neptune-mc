#ifndef NEPTUNE_BACKGROUND_HPP
#define NEPTUNE_BACKGROUND_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/vec3.hpp"
#include "neptune/physics/species.hpp"

namespace neptune {

struct BackgroundOneSpecies {
    vector<real_t> density;      // density of the species in each cell (m^-3)
    vector<real_t> velocity_x;   // velocity of the species in each cell in x
                                 // direction (m/s)
    vector<real_t> velocity_y;   // velocity of the species in each cell in y
                                 // direction (m/s)
    vector<real_t> velocity_z;   // velocity of the species in each cell in z
                                 // direction (m/s)
    vector<real_t> temperature;  // temperature of the species in each cell (eV)
};

/**
 * @brief This class contains the background plasma properties
 *      - density
 *      - velocity
 *      - temperature
 */
class Background {
   public:
    Background() = default;

    real_t density(Species species, index_t cell_index) const
    {
        return data_per_species[static_cast<index_t>(species)]
            .density[cell_index];
    }

    real_t velocity_x(Species species, index_t cell_index) const
    {
        return data_per_species[static_cast<index_t>(species)]
            .velocity_x[cell_index];
    }

    real_t velocity_y(Species species, index_t cell_index) const
    {
        return data_per_species[static_cast<index_t>(species)]
            .velocity_y[cell_index];
    }

    real_t velocity_z(Species species, index_t cell_index) const
    {
        return data_per_species[static_cast<index_t>(species)]
            .velocity_z[cell_index];
    }

    Vec3 velocity(Species species, index_t cell_index) const
    {
        return Vec3{
            data_per_species[static_cast<index_t>(species)]
                .velocity_x[cell_index],
            data_per_species[static_cast<index_t>(species)]
                .velocity_y[cell_index],
            data_per_species[static_cast<index_t>(species)]
                .velocity_z[cell_index],
        };
    }

    real_t temperature(Species species, index_t cell_index) const
    {
        return data_per_species[static_cast<index_t>(species)]
            .temperature[cell_index];
    }

    void set_species_data(Species species,
                          vector<real_t>&& density,
                          vector<real_t>&& velocity_x,
                          vector<real_t>&& velocity_y,
                          vector<real_t>&& velocity_z,
                          vector<real_t>&& temperature)
    {
        data_per_species[static_cast<index_t>(species)].density =
            std::move(density);
        data_per_species[static_cast<index_t>(species)].velocity_x =
            std::move(velocity_x);
        data_per_species[static_cast<index_t>(species)].velocity_y =
            std::move(velocity_y);
        data_per_species[static_cast<index_t>(species)].velocity_z =
            std::move(velocity_z);
        data_per_species[static_cast<index_t>(species)].temperature =
            std::move(temperature);
    }

   private:
    BackgroundOneSpecies data_per_species[NUM_SPECIES];
};

}  // namespace neptune

#endif