#include "neptune/estimation/estimatorgroup.hpp"

namespace neptune {

EstimatorGroup::EstimatorGroup(
    std::vector<std::tuple<Species, QOI, EstimatorType, bool>>
        estimators_definitions,
    index_t n_cells)
    : estimators()
{
    for (const auto& estimator : estimators_definitions) {
        estimators.emplace_back(std::get<0>(estimator), std::get<1>(estimator),
                                std::get<2>(estimator), n_cells,
                                std::get<3>(estimator));
    }

    // Fill the track-length estimator pointers
    for (auto& estimator : estimators) {
        if (estimator.estimator_type == EstimatorType::TrackLength) {
            track_length_estimators_by_species[static_cast<size_t>(
                                                   estimator.species)]
                .push_back(&estimator);
        }
        else if (estimator.estimator_type == EstimatorType::FinalTime) {
            final_estimators_by_species[static_cast<size_t>(estimator.species)]
                .push_back(&estimator);
        }
        else if (estimator.estimator_type == EstimatorType::KDMCEstimator) {
            kdmc_estimators_by_species[static_cast<size_t>(estimator.species)]
                .push_back(&estimator);
        }
    }
}

void EstimatorGroup::tally_track_length(
    const Particle& p,
    const Particle& p_old,
    const Mesh& mesh,
    const Background& background,
    const vector<std::unique_ptr<IonisationCollision>>& ionisation_collisions)
{
    if (track_length_estimators_by_species[static_cast<size_t>(p.species)]
            .empty()) {
        // No track-length estimators for this species, so we can skip the
        // calculations
        return;
    }

    // Calculate the total ionisation rate for this particle
    real_t R_i = 0;
    for (size_t i = 0; i < ionisation_collisions.size(); ++i) {
        real_t collision_rate =
            ionisation_collisions[i]->collision_rate(p_old, background);
        R_i += collision_rate;
    }

    real_t volume = mesh.get_volume(p_old.geometry.cell_id);
    real_t dt = p.t - p_old.t;
    real_t s;
    if (R_i * dt < 1.0e-4) {
        // Avoid division by zero if ionisation rate is very low
        s = p_old.weight * dt / volume * p_old.intensity;
    }
    else {
        s = p_old.weight * (1 - std::exp(-R_i * dt)) / R_i / volume *
            p_old.intensity;
    }

    // Loop over the relevant track-length estimators for this particle and
    // tally them
    for (auto* estimator :
         track_length_estimators_by_species[static_cast<size_t>(
             p_old.species)]) {
        estimator->tally(p_old, s, mesh.is_toroidal());
    }
}

void EstimatorGroup::tally_final(const Particle& p, const Mesh& mesh)
{
    if (final_estimators_by_species[static_cast<size_t>(p.species)].empty()) {
        // No final-time estimators for this species, so we can skip the
        // calculations
        return;
    }

    real_t volume = mesh.get_volume(p.geometry.cell_id);
    real_t weight = p.weight / volume * p.intensity;

    // Loop over the relevant final-time estimators for this particle and tally
    // them
    for (auto* estimator :
         final_estimators_by_species[static_cast<size_t>(p.species)]) {
        estimator->tally(p, weight, mesh.is_toroidal());
    }
}

void EstimatorGroup::tally_fluid(const Particle& p,
                                 const Particle& p_old,
                                 const Mesh& mesh)
{
    if (kdmc_estimators_by_species[static_cast<size_t>(p.species)].empty()) {
        // No fluid estimators for this species, so we can skip the calculations
        return;
    }

    real_t volume = mesh.get_volume(p.geometry.cell_id);
    real_t weight = p.weight / volume * p.intensity;

    real_t volume_old = mesh.get_volume(p_old.geometry.cell_id);
    real_t weight_old = p_old.weight / volume_old * p_old.intensity;

    // Loop over the relevant fluid estimators for this particle and tally them
    for (auto* estimator :
         kdmc_estimators_by_species[static_cast<size_t>(p.species)]) {
        if (estimator->qoi == QOI::KDMC_Fluid_Start_Density) {
            estimator->tally(p_old, weight_old, mesh.is_toroidal());
        }
        else if (estimator->qoi == QOI::KDMC_Fluid_End_Density) {
            estimator->tally(p, weight, mesh.is_toroidal());
        }
    }
}

void EstimatorGroup::combine(const EstimatorGroup& other)
{
    assert(estimators.size() == other.estimators.size());
    for (index_t i = 0; i < estimators.size(); ++i) {
        estimators[i].combine(other.estimators[i]);
    }
}

void EstimatorGroup::finalize_batch()
{
    for (auto& estimator : estimators) {
        estimator.finalize_batch();
    }
}

void EstimatorGroup::finalize()
{
    for (auto& estimator : estimators) {
        estimator.finalize();
    }
}

EstimatorGroup EstimatorGroup::make_new()
{
    std::vector<std::tuple<Species, QOI, EstimatorType, bool>>
        estimator_definitions;
    for (const auto& estimator : estimators) {
        estimator_definitions.emplace_back(estimator.species, estimator.qoi,
                                           estimator.estimator_type,
                                           estimator.estimate_variance);
    }
    return EstimatorGroup(estimator_definitions,
                          estimators[0].tally_vector.size());
}

}  // namespace neptune