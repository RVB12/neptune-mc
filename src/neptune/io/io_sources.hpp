#ifndef NEPTUNE_IO_SOURCES_HPP
#define NEPTUNE_IO_SOURCES_HPP

#include "neptune/io/io_common.hpp"
#include "neptune/sources/recombinationsource.hpp"
#include "neptune/sources/source.hpp"
#include "neptune/sources/volumesource.hpp"
#include "neptune/sources/wallreflectionsource.hpp"

namespace neptune {

inline vector<std::unique_ptr<Source>> read_sources(const H5::Group& sources_group,
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
    std::string stratification_method =
        read_string(sources_group.openDataSet("stratification"));
    if (stratification_method == "source_strength") {
        real_t total_intensity = 0;
        for (const auto& source : sources) {
            total_intensity += source->intensity();
        }
        for (auto& source : sources) {
            source->num_particles_ =
                std::ceil(n_particles * source->intensity() / total_intensity);
        }
    }
    else if (stratification_method == "equal") {
        for (auto& source : sources) {
            source->num_particles_ = static_cast<index_t>(
                n_particles /
                sources.size());  // equal number of particles for each source
        }
    }
    else {
        throw std::runtime_error("Unknown stratification method: " +
                                 stratification_method);
    }

    // Convert std::vector to neptune::vector
    vector<std::unique_ptr<Source>> sources_vector(sources.size());
    for (size_t i = 0; i < sources.size(); i++) {
        sources_vector[i] = std::move(sources[i]);
    }
    return sources_vector;
}
}  // namespace neptune

#endif