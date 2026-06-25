#include "neptune/boundaryconditions/trim.hpp"

#include <fstream>
#include <stdexcept>

namespace neptune {

TrimReflectionData::TrimReflectionData(
    Species incident_species,
    vector<real_t>& incident_energy,
    vector<real_t>& incident_polar_angle,
    vector<real_t>& reflection_coefficient,
    vector<real_t>& reflected_energy,
    vector<real_t>& reflected_polar_angle,
    vector<real_t>& reflected_azimuthal_angle)
    : indicent_species(incident_species),
      incident_energy(incident_energy),
      incident_polar_angle(incident_polar_angle),
      reflection_coefficient(reflection_coefficient),
      reflected_energy(reflected_energy),
      reflected_polar_angle(reflected_polar_angle),
      reflected_azimuthal_angle(reflected_azimuthal_angle)
{
    // Check that the dimensions of the data are consistent
    index_t n_tables = incident_energy.size();
    if (incident_polar_angle.size() != n_tables ||
        reflection_coefficient.size() != n_tables) {
        throw std::runtime_error(
            "Incident energy, polar angle and reflection coefficient "
            "must have the same size");
    }
    index_t n_energy_bins = reflected_energy.size() / n_tables;
    if (reflected_energy.size() != n_tables * n_energy_bins) {
        throw std::runtime_error(
            "Reflected energy must have size n_tables x n_energy_bins");
    }
    index_t n_polar_bins =
        reflected_polar_angle.size() / (n_tables * n_energy_bins);
    if (reflected_polar_angle.size() !=
        n_tables * n_energy_bins * n_polar_bins) {
        throw std::runtime_error(
            "Reflected polar angle must have size "
            "n_tables x n_energy_bins x n_polar_bins");
    }
    index_t n_azimuthal_bins =
        reflected_azimuthal_angle.size() / (n_tables * n_energy_bins);
    if (reflected_azimuthal_angle.size() !=
        n_tables * n_energy_bins * n_azimuthal_bins) {
        throw std::runtime_error(
            "Reflected azimuthal angle must have size "
            "n_tables x n_energy_bins x n_azimuthal_bins");
    }
}

TrimReflectionData::TrimReflectionData(
    Species incident_species,
    vector<real_t>&& incident_energy,
    vector<real_t>&& incident_polar_angle,
    vector<real_t>&& reflection_coefficient,
    vector<real_t>&& reflected_energy,
    vector<real_t>&& reflected_polar_angle,
    vector<real_t>&& reflected_azimuthal_angle)
    : indicent_species(incident_species),
      incident_energy(std::move(incident_energy)),
      incident_polar_angle(std::move(incident_polar_angle)),
      reflection_coefficient(std::move(reflection_coefficient)),
      reflected_energy(std::move(reflected_energy)),
      reflected_polar_angle(std::move(reflected_polar_angle)),
      reflected_azimuthal_angle(std::move(reflected_azimuthal_angle))
{
    // Check that the dimensions of the data are consistent
    index_t n_tables = incident_energy.size();
    if (incident_polar_angle.size() != n_tables ||
        reflection_coefficient.size() != n_tables) {
        throw std::runtime_error(
            "Incident energy, polar angle and reflection coefficient "
            "must have the same size");
    }
    index_t n_energy_bins = reflected_energy.size() / n_tables;
    if (reflected_energy.size() != n_tables * n_energy_bins) {
        throw std::runtime_error(
            "Reflected energy must have size n_tables x n_energy_bins");
    }
    index_t n_polar_bins =
        reflected_polar_angle.size() / (n_tables * n_energy_bins);
    if (reflected_polar_angle.size() !=
        n_tables * n_energy_bins * n_polar_bins) {
        throw std::runtime_error(
            "Reflected polar angle must have size "
            "n_tables x n_energy_bins x n_polar_bins");
    }
    index_t n_azimuthal_bins =
        reflected_azimuthal_angle.size() / (n_tables * n_energy_bins);
    if (reflected_azimuthal_angle.size() !=
        n_tables * n_energy_bins * n_azimuthal_bins) {
        throw std::runtime_error(
            "Reflected azimuthal angle must have size "
            "n_tables x n_energy_bins x n_azimuthal_bins");
    }
}

Vec3 TrimReflectionData::sample_reflection(Particle& p,
                                           Vec3 surface_normal) const
{
    // Calculate the incident energy and polar angle of the particle
    real_t energy = 0.5 * p.geometry.velocity.norm2() *
                    SPECIES_MASS[static_cast<size_t>(p.species)] *
                    constants::m_u / constants::eV_to_Joules;

    real_t polar_angle = std::acos(-p.geometry.velocity.dot(surface_normal) /
                                   p.geometry.velocity.norm());
    polar_angle = polar_angle * 180 / M_PI;

    // Find the index in the incident energy and polar angle arrays
    index_t energy_index = binary_search(incident_energy, energy);
    index_t polar_index = binary_search(incident_polar_angle, polar_angle);

    // Now sample three uniformly distributed random numbers to determine
    // which bin of the reflected energy, polar angle and azimuthal angle to
    // use
    std::uniform_real_distribution<real_t> dist_u(0, 1);

    real_t u_energy = dist_u(p.rng);
    real_t u_polar = dist_u(p.rng);
    real_t u_azimuthal = dist_u(p.rng);

    index_t n_tables = incident_energy.size();
    index_t n_energy_bins = reflected_energy.size() / n_tables;
    index_t n_polar_bins =
        reflected_polar_angle.size() / (n_tables * n_energy_bins);
    index_t n_azimuthal_bins =
        reflected_azimuthal_angle.size() / (n_tables * n_energy_bins);

    index_t reflected_energy_index = std::floor(u_energy * n_energy_bins);
    index_t reflected_polar_index = std::floor(u_polar * n_polar_bins);
    index_t reflected_azimuthal_index =
        std::floor(u_azimuthal * n_azimuthal_bins);

    // Get the reflected energy, polar angle and azimuthal angle from the
    // data
    // TODO: instead of just using the corresponding bin, we should
    // do multilinear interpolation between the bins to get a more accurate
    // result
    real_t reflected_energy_value =
        reflected_energy[energy_index * n_energy_bins + reflected_energy_index];
    real_t reflected_polar_cos =
        reflected_polar_angle[energy_index * n_energy_bins * n_polar_bins +
                              reflected_energy_index * n_polar_bins +
                              reflected_polar_index];
    real_t reflected_azimuthal_cos =
        reflected_azimuthal_angle[energy_index * n_energy_bins *
                                      n_azimuthal_bins +
                                  reflected_energy_index * n_azimuthal_bins +
                                  reflected_azimuthal_index];

    // Convert the reflected energy, polar angle and azimuthal angle to a
    // velocity vector
    real_t reflected_speed = std::sqrt(
        2 * reflected_energy_value * constants::eV_to_Joules /
        (SPECIES_MASS[static_cast<size_t>(p.species)] * constants::m_u));

    real_t reflected_polar_sin =
        std::sqrt(1 - reflected_polar_cos * reflected_polar_cos);
    real_t reflected_azimuthal_sin =
        std::sqrt(1 - reflected_azimuthal_cos * reflected_azimuthal_cos);

    Vec3 normal2 = p.geometry.velocity -
                   p.geometry.velocity.dot(surface_normal) * surface_normal;
    if (normal2.norm() > 1e-6) {
        normal2 /= normal2.norm();
    }
    else {
        // If the normal2 vector is very small, we can choose any vector
        // perpendicular to the surface normal
        if (std::abs(surface_normal.z) < 0.9) {
            normal2 = Vec3(0, 0, 1).cross(surface_normal);
        }
        else {
            normal2 = Vec3(1, 0, 0).cross(surface_normal);
        }
        normal2 /= normal2.norm();
    }
    Vec3 normal3 = surface_normal.cross(normal2);

    Vec3 reflected_velocity =
        reflected_speed *
        (surface_normal * reflected_polar_cos +
         normal2 * reflected_polar_sin * reflected_azimuthal_cos +
         normal3 * reflected_polar_sin * reflected_azimuthal_sin);

    return reflected_velocity;
};

TrimReflectionData read_trim_data(const std::string& filename)
{
    std::ifstream trim_file(filename);

    if (!trim_file.is_open()) {
        throw std::runtime_error("Could not open TRIM data file: " + filename);
    }

    // Ignore first 6 lines of the file, they contain metadata that we don't
    // need
    std::string line;
    for (int i = 0; i < 6; i++) {
        std::getline(trim_file, line);
    }

    // Hardcode the dimensions of the data for now
    // Should perhaps be read from the file in the future
    index_t n_tables = 12 * 7;
    index_t n_energy_bins = 5;
    index_t n_polar_bins = 5;
    index_t n_azimuthal_bins = 5;

    vector<real_t> incident_energy(n_tables);
    vector<real_t> incident_polar_angle(n_tables);
    vector<real_t> reflection_coefficient(n_tables);
    vector<real_t> reflected_energy(n_tables * n_energy_bins);
    vector<real_t> reflected_polar_angle(n_tables * n_energy_bins *
                                         n_polar_bins);
    vector<real_t> reflected_azimuthal_angle(n_tables * n_energy_bins * n_polar_bins *
                                             n_azimuthal_bins);

    for (index_t i = 0; i < n_tables; i++) {
        real_t incident_atomic_number;
        trim_file >> incident_atomic_number;
        real_t incident_atomic_mass;
        trim_file >> incident_atomic_mass;
        real_t wall_atomic_number;
        trim_file >> wall_atomic_number;
        real_t wall_atomic_mass;
        trim_file >> wall_atomic_mass;

        trim_file >> incident_energy[i];
        trim_file >> incident_polar_angle[i];
        trim_file >> reflection_coefficient[i];

        // skip the rest of the line
        std::getline(trim_file, line);
        // skip the next line too
        std::getline(trim_file, line);

        // Now read the reflected energy
        for (index_t j = 0; j < n_energy_bins; j++) {
            trim_file >> reflected_energy[i * n_energy_bins + j];
        }

        // skip the rest of the line
        std::getline(trim_file, line);
        // skip the next line too
        std::getline(trim_file, line);

        // Now read the reflected polar angle
        for (index_t j = 0; j < n_energy_bins; j++) {
            for (index_t jj = 0; jj < n_polar_bins; jj++) {
                trim_file >>
                    reflected_polar_angle[i * n_energy_bins * n_polar_bins +
                                          j * n_polar_bins + jj];
            }
            // skip the rest of the line
            std::getline(trim_file, line);
        }

        // skip the next line
        std::getline(trim_file, line);

        // Now read the reflected azimuthal angle
        for (index_t j = 0; j < n_energy_bins; j++) {
            for (index_t jj = 0; jj < n_polar_bins; jj++) {
                for (index_t jjj = 0; jjj < n_azimuthal_bins; jjj++) {
                    trim_file >>
                        reflected_azimuthal_angle[i * n_energy_bins *
                                                      n_polar_bins *
                                                      n_azimuthal_bins +
                                                  j * n_polar_bins *
                                                      n_azimuthal_bins +
                                                  jj * n_azimuthal_bins + jjj];
                }
                // skip the rest of the line
                std::getline(trim_file, line);
            }
        }
    }

    // TODO: determine species from the incident atomic number and mass
    Species incident_species = Species::D;

    return TrimReflectionData(
        incident_species, std::move(incident_energy),
        std::move(incident_polar_angle), std::move(reflection_coefficient),
        std::move(reflected_energy), std::move(reflected_polar_angle),
        std::move(reflected_azimuthal_angle));
}

}  // namespace neptune