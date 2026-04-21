#ifndef NEPTUNE_ESTIMATOR_HPP
#define NEPTUNE_ESTIMATOR_HPP

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/physics/background.hpp"
#include "neptune/physics/particle.hpp"

namespace neptune {

enum class QOI {
    Density,
    Momentum_x,
    Momentum_y,
    Momentum_z,
    Energy,
    Energy_x,
    Energy_y,
    Energy_z,
};

inline QOI stringtoQOI(std::string s)
{
    if (s == "density") return QOI::Density;
    if (s == "momentum_x") return QOI::Momentum_x;
    if (s == "momentum_y") return QOI::Momentum_y;
    if (s == "momentum_z") return QOI::Momentum_z;
    if (s == "energy") return QOI::Energy;
    if (s == "energy_x") return QOI::Energy_x;
    if (s == "energy_y") return QOI::Energy_y;
    if (s == "energy_z") return QOI::Energy_z;
    throw std::runtime_error("Unknown QOI: " + s);
}

enum class EstimatorType { TrackLength, FinalTime };

inline EstimatorType stringtoEstimatorType(std::string s)
{
    if (s == "track_length") return EstimatorType::TrackLength;
    if (s == "final_time") return EstimatorType::FinalTime;
    throw std::runtime_error("Unknown estimator type: " + s);
}

/**
 * An estimator for a quantity of interest (QOI).
 */
class Estimator {
   public:
    Estimator(Species species,
              QOI qoi,
              EstimatorType estimator_type,
              index_t n_cells,
              bool estimate_variance);

    /**
     * @brief Tally the contribution
     *
     * @param p The particle to tally.
     *          Note, for track-length estimators, this should be the particle
     * before the move For final-time estimators, this should be the particle at
     * the final time of the simulation
     * @param weight The weight to use for this contribution
     *               For Track-length, this should be \int_t^t+time_in_cell
     * weight(t) dt * intensity / volume For Final-time, this should be
     * weight(t_final) * intensity / volume
     */
    void tally(const Particle& p, real_t weighted);

    void combine(const Estimator& other);

    const bool estimate_variance;
    const Species species;
    const QOI qoi;
    const EstimatorType estimator_type;
    vector<real_t> tally_vector;
    vector<real_t> tally_mean;
    vector<real_t> tally_variance;
};

/**
 * A group of estimators for different quantities of interest (QOI) and/or
 * different types of estimators (e.g. track-length and final-time).
 *
 * We specifically use a class for this so that some computations can be shared
 * between the different estimators, e.g. the computation of the weight
 * reduction for ionisation.
 */
class EstimatorGroup {
   public:
    EstimatorGroup(std::vector<std::tuple<Species, QOI, EstimatorType>> estimators_definitions, index_t n_cells);

    void tally_track_length(const Particle& p,
                            const Particle& p_old,
                            const Mesh& mesh,
                            const Background& background,
                            const vector<std::unique_ptr<IonisationCollision>>&
                                ionisation_collisions);

    void tally_final(const Particle& p, const Mesh& mesh);

    // Combine the estimates of two estimator groups
    void combine(const EstimatorGroup& other);

    // Make a new estimator group with the same estimators but all estimates
    // reset to zero used for multi-threading
    EstimatorGroup make_new();

   private:
    std::vector<Estimator> estimators;

    // An array of pointers to the track-length estimators, indexed by species, used
    // to quickly find the relevant estimators for a given particle
    std::vector<Estimator*> track_length_estimators_by_species[NUM_SPECIES];
    std::vector<Estimator*> final_estimators_by_species[NUM_SPECIES];
};

#pragma omp declare reduction(EstimatorGroupReduction \
:EstimatorGroup : omp_out.combine(omp_in))            \
    initializer(omp_priv = omp_orig.make_new())

}  // namespace neptune

#endif  // NEPTUNE_ESTIMATOR_HPP