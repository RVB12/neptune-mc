#ifndef NEPTUNE_IO_BOUNDARYCONDITIONS_HPP
#define NEPTUNE_IO_BOUNDARYCONDITIONS_HPP

#include "neptune/io/io_common.hpp"
#include "neptune/boundaryconditions/boundarycondition.hpp"
#include "neptune/boundaryconditions/trim.hpp"

namespace neptune {

inline vector<std::unique_ptr<BoundaryCondition>> read_boundary_conditions(
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
        else if( bc_type == "diffuse") {
            real_t fixed_energy = read_real_number(
                bc_group.openDataSet("fixed_energy"));
            real_t absorption_coefficient = read_real_number(
                bc_group.openDataSet("absorption_coefficient"));
            boundary_conditions.push_back(
                std::make_unique<DiffuseBC>(fixed_energy, absorption_coefficient));
        }
        else if (bc_type == "trim") {
            std::string trim_filename =
                read_string(bc_group.openDataSet("trim_filename"));
            boundary_conditions.push_back(
                std::make_unique<TrimReflection>(trim_filename));
        }
        else if (bc_type == "transparent") {
            boundary_conditions.push_back(std::make_unique<TransparentBC>());
        }
        else if (bc_type == "periodic") {
            throw std::runtime_error("Periodic boundary conditions are not supported yet");
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
}  // namespace neptune

#endif