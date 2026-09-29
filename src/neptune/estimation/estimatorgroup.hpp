#ifndef NEPTUNE_ESTIMATORGROUP_HPP
#define NEPTUNE_ESTIMATORGROUP_HPP

#include <vector>

#include "neptune/estimation/estimator.hpp"

namespace neptune {

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
    EstimatorGroup(std::vector<std::tuple<Species, QOI, EstimatorType, bool, bool, index_t>>
                       estimators_definitions,
                   index_t n_cells);

    void tally_track_length(const Particle& p,
                            const Particle& p_old,
                            const Mesh& mesh,
                            const Background& background,
                            const vector<std::unique_ptr<IonisationCollision>>&
                                ionisation_collisions);

    void tally_final(const Particle& p, const Mesh& mesh);

    void tally_fluid(const Particle& p, const Particle& p_old, const Mesh& mesh);

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
    std::vector<Estimator*> kdmc_estimators_by_species[NUM_SPECIES];
};

#pragma omp declare reduction(                                            \
        EstimatorGroupReduction:EstimatorGroup : omp_out.combine(omp_in)) \
    initializer(omp_priv = omp_orig.make_new())

}  // namespace neptune

#endif