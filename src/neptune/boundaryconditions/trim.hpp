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

    std::pair<real_t, Vec3> sample_reflection(Particle& p,
                                              Vec3 surface_normal) const;

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
    TrimReflection(const std::string& filename,
                   real_t E_wall,
                   real_t absorption_coefficient,
                   real_t fast_reflection_multiplier)
        : data(read_trim_data(filename)),
          E_wall(E_wall),
          absorption_coefficient(absorption_coefficient),
          fast_reflection_multiplier(fast_reflection_multiplier)
    {}
    ~TrimReflection() = default;

    void apply(Particle& p,
               const Mesh& mesh,
               const CellCrossing& crossing) override
    {
        Vec3 normal = mesh.get_surface_normal(
            p.geometry.cell_id, p.geometry.boundary_id, p.geometry.position);
        // Determine if the particle undergoes fast reflection or thermal
        // reflection
        auto [probability, velocity] = data.sample_reflection(p, normal);

        real_t Rt = 1.0 - absorption_coefficient;
        real_t Rfast = std::min(fast_reflection_multiplier * probability, Rt);
        real_t Rthermal = Rt - Rfast;
        // Determine if we should keep the sampled velocity from trim, or if we
        // should sample a thermal velocity

        std::uniform_real_distribution<real_t> dist_u(0, 1);
        real_t u = dist_u(p.rng);
        if (u * Rt < Rthermal) {
            // Sample thermal velocity from cosine distribution
            Vec3 normal = mesh.get_surface_normal(p.geometry.cell_id,
                                                  p.geometry.boundary_id,
                                                  p.geometry.position);
            // Generate the direction of the reflected particle using a cosine
            // distribution
            std::uniform_real_distribution<real_t> dist_u(0, 1);
            real_t u = dist_u(p.rng);
            Vec3 tangent1;
            if (std::abs(normal.z) < 0.9) {
                tangent1 = Vec3(0, 0, 1).cross(normal);
            }
            else {
                tangent1 = Vec3(1, 0, 0).cross(normal);
            }
            tangent1 /= tangent1.norm();
            Vec3 tangent2 = normal.cross(tangent1);

            real_t phi = 2 * M_PI * dist_u(p.rng);
            real_t cos_theta = std::sqrt(dist_u(p.rng));
            real_t sin_theta = std::sqrt(1 - cos_theta * cos_theta);

            Vec3 direction = cos_theta * normal +
                             sin_theta * std::cos(phi) * tangent1 +
                             sin_theta * std::sin(phi) * tangent2;
            direction /= direction.norm();

            // Set the speed of the particle to the fixed energy
            // Note, I was initially confused about the eirene documentation
            // It says EWALL = 𝐸𝑚𝑒𝑎𝑛 = 1.5 · 𝑇𝑊
            // but the 1.5 factor is already included in the definition of the temperature
            real_t fixed_energy = E_wall;
            real_t mass =
                constants::m_u * SPECIES_MASS[static_cast<size_t>(p.species)];
            real_t velocity_magnitude =
                std::sqrt(2 * fixed_energy * constants::eV_to_Joules / mass);

            p.geometry.velocity = direction * velocity_magnitude;
        }
        else {
            p.geometry.velocity = velocity;
        }

        p.weight *= Rt;
    }

   private:
    const TrimReflectionData data;
    const real_t E_wall;
    const real_t fast_reflection_multiplier;
    const real_t absorption_coefficient;
};

}  // namespace neptune

#endif  // NEPTUNE_TRIM_HPP