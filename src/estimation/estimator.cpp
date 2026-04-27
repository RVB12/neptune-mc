#include "neptune/estimation/estimator.hpp"

#include "neptune/geometry/cart_to_tor.hpp"
#include "neptune/physics/physicalconstants.hpp"
namespace neptune {

Estimator::Estimator(Species species,
                     QOI qoi,
                     EstimatorType estimator_type,
                     index_t n_cells,
                     bool estimate_variance)
    : species(species),
      qoi(qoi),
      estimator_type(estimator_type),
      estimate_variance(estimate_variance),
      tally_vector(n_cells),
      tally_mean(estimate_variance ? n_cells : 0),
      tally_variance(estimate_variance ? n_cells : 0)
{
    for (index_t i = 0; i < n_cells; ++i) {
        tally_vector[i] = 0;
    }
    if (estimate_variance) {
        for (index_t i = 0; i < n_cells; ++i) {
            tally_mean[i] = 0;
            tally_variance[i] = 0;
        }
    }
}

void Estimator::combine(const Estimator& other)
{
    assert(qoi == other.qoi);
    assert(estimator_type == other.estimator_type);
    assert(tally_vector.size() == other.tally_vector.size());
    // Note, in principle, the number of samples should really be the same
    // but this may not be the case if the counter has been reset
    // in the reduction
    n_samples = std::max(n_samples, other.n_samples);
    for (index_t i = 0; i < tally_vector.size(); ++i) {
        tally_vector[i] += other.tally_vector[i];
    }
    if (estimate_variance) {
        for (index_t i = 0; i < tally_vector.size(); ++i) {
            tally_mean[i] += other.tally_mean[i];
            tally_variance[i] += other.tally_variance[i];
        }
    }
}

void Estimator::tally(const Particle& p, real_t weight, bool is_toroidal)
{
    switch (qoi) {
        case QOI::Density:
            tally_vector[p.geometry.cell_id] += weight;
            break;
        case QOI::Velocity_x:
            if (is_toroidal) {
                // Note, this calculation is not perfect,
                // even if the velocity is constant, the toroidal velocity
                // can still vary so we should actually integrate
                Vec3 vel_tor = transform_cart_to_tor(p.geometry.position,
                                                     p.geometry.velocity);
                tally_vector[p.geometry.cell_id] += vel_tor.x * weight;
            }
            else {
                tally_vector[p.geometry.cell_id] +=
                    p.geometry.velocity.x * weight;
            }
            break;
        case QOI::Velocity_y:
            if (is_toroidal) {
                // Note, this calculation is not perfect,
                // even if the velocity is constant, the toroidal velocity
                // can still vary so we should actually integrate
                Vec3 vel_tor = transform_cart_to_tor(p.geometry.position,
                                                     p.geometry.velocity);
                tally_vector[p.geometry.cell_id] += vel_tor.y * weight;
            }
            else {
                tally_vector[p.geometry.cell_id] +=
                    p.geometry.velocity.y * weight;
            }
            break;
        case QOI::Velocity_z:
            tally_vector[p.geometry.cell_id] += p.geometry.velocity.z * weight;
            break;
        case QOI::Energy:
            tally_vector[p.geometry.cell_id] +=
                (0.5 * p.geometry.velocity.norm2()) * weight;
            break;
        case QOI::Energy_x:
            if (is_toroidal) {
                // Note, this calculation is not perfect,
                // even if the velocity is constant, the toroidal velocity
                // can still vary so we should actually integrate
                Vec3 vel_tor = transform_cart_to_tor(p.geometry.position,
                                                     p.geometry.velocity);
                tally_vector[p.geometry.cell_id] +=
                    (0.5 * vel_tor.x * vel_tor.x) *
                    weight;
            }
            else {
                tally_vector[p.geometry.cell_id] +=
                    (0.5 * p.geometry.velocity.x * p.geometry.velocity.x) *
                    weight;
            }
            break;
        case QOI::Energy_y:
            if (is_toroidal) {
                // Note, this calculation is not perfect,
                // even if the velocity is constant, the toroidal velocity
                // can still vary so we should actually integrate
                Vec3 vel_tor = transform_cart_to_tor(p.geometry.position,
                                                     p.geometry.velocity);
                tally_vector[p.geometry.cell_id] +=
                    (0.5 * vel_tor.y * vel_tor.y) *
                    weight;
            }
            else {
                tally_vector[p.geometry.cell_id] +=
                    (0.5 * p.geometry.velocity.y * p.geometry.velocity.y) *
                    weight;
            }
            break;
        case QOI::Energy_z:
            tally_vector[p.geometry.cell_id] +=
                (0.5 * p.geometry.velocity.z * p.geometry.velocity.z) * weight;
            break;
    }
}

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
        real_t collision_rate = ionisation_collisions[i]->collision_rate(
            background.density(ionisation_collisions[i]->background_species(),
                               p_old.geometry.cell_id),
            background.velocity(ionisation_collisions[i]->background_species(),
                                p_old.geometry.cell_id),
            background.temperature(
                ionisation_collisions[i]->background_species(),
                p_old.geometry.cell_id),
            p_old);
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