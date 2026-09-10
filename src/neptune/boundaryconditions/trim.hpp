#ifndef NEPTUNE_TRIM_HPP
#define NEPTUNE_TRIM_HPP

#include "neptune/boundaryconditions/boundarycondition.hpp"
#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/physics/physicalconstants.hpp"
#include "neptune/physics/species.hpp"

namespace neptune {

// This class stores the data for the TRIM reflection model
//
class TrimReflectionData {
   public:
    TrimReflectionData(Species incident_species,
                       vector<real_t>& incident_energy,
                       vector<real_t>& incident_polar_angle,
                       vector<real_t>& reflection_coefficient,
                       vector<real_t>& reflected_energy,
                       vector<real_t>& reflected_polar_angle,
                       vector<real_t>& reflected_azimuthal_angle);

    TrimReflectionData(Species incident_species,
                       vector<real_t>&& incident_energy,
                       vector<real_t>&& incident_polar_angle,
                       vector<real_t>&& reflection_coefficient,
                       vector<real_t>&& reflected_energy,
                       vector<real_t>&& reflected_polar_angle,
                       vector<real_t>&& reflected_azimuthal_angle);

    std::pair<real_t, Vec3> sample_reflection(Particle& p, Vec3 surface_normal) const;

   private:
    // The species of the incident particle
    Species indicent_species;
    // incident energy of the particle (eV)
    // dimension: n_tables
    vector<real_t> incident_energy;
    // polar angle of the incident particle (degrees)
    // dimension: n_tables
    vector<real_t> incident_polar_angle;
    // Reflection coefficient, i.e. the fraction of particles reflected
    // dimension: n_tables
    vector<real_t> reflection_coefficient;
    // Energy of the reflected particle (eV)
    // dimension: n_tables x n_energy_bins
    vector<real_t> reflected_energy;
    // Cosine of the polar angle of the reflected particle
    // dimension: n_tables x n_energy_bins x n_polar_bins
    vector<real_t> reflected_polar_angle;
    // Cosine of the azimuthal angle of the reflected particle
    // dimension: n_tables x n_energy_bins x n_azimuthal_bins
    vector<real_t> reflected_azimuthal_angle;
};

TrimReflectionData read_trim_data(const std::string& filename);

class TrimReflection : public BoundaryCondition {
   public:
    TrimReflection(const std::string& filename) : data(read_trim_data(filename)) {}
    ~TrimReflection() = default;

    void apply(Particle& p, const Mesh& mesh, const CellCrossing& crossing) override
    {
        Vec3 normal = mesh.get_surface_normal(p.geometry.cell_id, p.geometry.boundary_id, p.geometry.position);
        auto [probability, velocity] = data.sample_reflection(p, normal);
        p.geometry.velocity = velocity;
        p.weight *= probability;
    }

   private:
    const TrimReflectionData data;
};

}  // namespace neptune

#endif  // NEPTUNE_TRIM_HPP