#ifndef NEPTUNE_ESTIMATOR_HPP
#define NEPTUNE_ESTIMATOR_HPP

#include "neptune/collisions/abstract_collisions.hpp"
#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/estimation/qoi.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/physics/background.hpp"
#include "neptune/physics/particle.hpp"

namespace neptune {

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
     * @param is_toroidal Whether the velocity should be transformed from
     * Carthesian to toroidal coordinates
     */
    void tally(const Particle& p, real_t weighted, bool is_toroidal);

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

}  // namespace neptune

#endif  // NEPTUNE_ESTIMATOR_HPP