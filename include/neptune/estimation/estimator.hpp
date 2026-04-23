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

inline std::string QOItostring(QOI qoi)
{
    switch (qoi) {
        case QOI::Density:
            return "density";
        case QOI::Momentum_x:
            return "momentum_x";
        case QOI::Momentum_y:
            return "momentum_y";
        case QOI::Momentum_z:
            return "momentum_z";
        case QOI::Energy:
            return "energy";
        case QOI::Energy_x:
            return "energy_x";
        case QOI::Energy_y:
            return "energy_y";
        case QOI::Energy_z:
            return "energy_z";
        default:
            throw std::runtime_error("Unknown QOI");
    }
}

enum class EstimatorType { TrackLength, FinalTime };

inline EstimatorType stringtoEstimatorType(std::string s)
{
    if (s == "track_length") return EstimatorType::TrackLength;
    if (s == "final_time") return EstimatorType::FinalTime;
    throw std::runtime_error("Unknown estimator type: " + s);
}

inline std::string estimator_type_tostring(EstimatorType estimator_type)
{
    switch (estimator_type) {
        case EstimatorType::TrackLength:
            return "track_length";
        case EstimatorType::FinalTime:
            return "final_time";
        default:
            throw std::runtime_error("Unknown estimator type");
    }
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

    vector<real_t>& get_estimate()
    {
        if (estimate_variance) {
            return tally_mean;
        }
        else {
            return tally_vector;
        }
    }

    vector<real_t>& get_variance()
    {
        if (estimate_variance) {
            return tally_variance;
        }
        else {
            throw std::runtime_error(
                "Variance estimation is not enabled for this estimator");
        }
    }

    /**
     * After running a batch, finalize the estimates
     * For the normal tally, nothing needs to be done
     * For the variance estimation, we update the mean and variance
     */
    void finalize_batch()
    {
        n_samples += 1;
        if (estimate_variance) {
            for (index_t i = 0; i < tally_vector.size(); i++) {
                real_t mean = tally_mean[i];
                real_t M2 = tally_variance[i];
                real_t new_value = tally_vector[i];
                real_t delta = new_value - mean;
                mean += delta / n_samples;
                M2 += delta * (new_value - mean);
                tally_mean[i] = mean;
                tally_variance[i] = M2;
                tally_vector[i] =
                    0;  // reset the tally vector for the next batch
            }
        }
    }

    void finalize()
    {
        if (estimate_variance) {
            for (index_t i = 0; i < tally_variance.size(); i++) {
                tally_variance[i] /= n_samples - 1;
            }
        }
        else {
            if (n_samples > 1) {
                for (index_t i = 0; i < tally_vector.size(); i++) {
                    tally_vector[i] /= n_samples;
                }
            }
        }
    }

    const bool estimate_variance;
    const Species species;
    const QOI qoi;
    const EstimatorType estimator_type;
    vector<real_t> tally_vector;
    vector<real_t> tally_mean;
    vector<real_t> tally_variance;
    index_t n_samples = 0;
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
    EstimatorGroup(std::vector<std::tuple<Species, QOI, EstimatorType, bool>>
                       estimators_definitions,
                   index_t n_cells);

    void tally_track_length(const Particle& p,
                            const Particle& p_old,
                            const Mesh& mesh,
                            const Background& background,
                            const vector<std::unique_ptr<IonisationCollision>>&
                                ionisation_collisions);

    void tally_final(const Particle& p, const Mesh& mesh);

    // Combine the estimates of two estimator groups
    void combine(const EstimatorGroup& other);

    void finalize_batch();

    void finalize();

    // Make a new estimator group with the same estimators but all estimates
    // reset to zero used for multi-threading
    EstimatorGroup make_new();

    std::vector<Estimator> estimators;

    // An array of pointers to the track-length estimators, indexed by species,
    // used to quickly find the relevant estimators for a given particle
    std::vector<Estimator*> track_length_estimators_by_species[NUM_SPECIES];
    std::vector<Estimator*> final_estimators_by_species[NUM_SPECIES];
};

#pragma omp declare reduction(                                            \
        EstimatorGroupReduction:EstimatorGroup : omp_out.combine(omp_in)) \
    initializer(omp_priv = omp_orig.make_new())

}  // namespace neptune

#endif  // NEPTUNE_ESTIMATOR_HPP