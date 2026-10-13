#include "neptune/estimation/estimator.hpp"

#include "neptune/geometry/cart_to_tor.hpp"
#include "neptune/physics/physicalconstants.hpp"
namespace neptune {

Estimator::Estimator(Species species,
                     QOI qoi,
                     EstimatorType estimator_type,
                     index_t n_cells,
                     bool estimate_variance,
                     bool keep_all_batches,
                     index_t n_batches)
    : species(species),
      qoi(qoi),
      estimator_type(estimator_type),
      estimate_variance(estimate_variance),
      keep_all_batches(keep_all_batches),
      tally_vector(n_cells),
      tally_mean(estimate_variance ? n_cells : 0),
      tally_variance(estimate_variance ? n_cells : 0),
      all_batches(keep_all_batches ? n_cells * n_batches : 0),
      n_batches(n_batches),
      KDMC_avg_fluid_step_duration_all_batches(keep_all_batches ? n_batches : 0)
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
    if (keep_all_batches) {
        // Initializing is not strictly necessary, since it will be overwritten
        // later
        for (index_t i = 0; i < n_cells * n_batches; ++i) {
            all_batches[i] = 0;
        }
        for (index_t i = 0; i < n_batches; ++i) {
            KDMC_avg_fluid_step_duration_all_batches[i] = 0;
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
    i_batch = std::max(i_batch, other.i_batch);
    for (index_t i = 0; i < tally_vector.size(); ++i) {
        tally_vector[i] += other.tally_vector[i];
    }
    if (estimate_variance) {
        for (index_t i = 0; i < tally_vector.size(); ++i) {
            tally_mean[i] += other.tally_mean[i];
            tally_variance[i] += other.tally_variance[i];
        }
    }
    if (keep_all_batches) {
        assert(all_batches.size() == other.all_batches.size());
        for (index_t i = 0; i < all_batches.size(); ++i) {
            all_batches[i] += other.all_batches[i];
        }
    }
    if(qoi == QOI::KDMC_Fluid_Start_Density || qoi == QOI::KDMC_Fluid_End_Density) {
        KDMC_avg_fluid_step_duration += other.KDMC_avg_fluid_step_duration;
        KDMC_nb_fluid_steps += other.KDMC_nb_fluid_steps;
    }
}

void Estimator::tally(const Particle& p,
                      const Particle& p_new,
                      real_t weight,
                      bool is_toroidal)
{
    switch (qoi) {
        case QOI::Density:
            tally_vector[p.geometry.cell_id] += weight;
            break;
        case QOI::Velocity_x:
            if (is_toroidal) {
                if (p_new.t != p.t) {
                    // The integral int_0^s v_r dt, with v_r the radial
                    // velocity, is equal to: int_0^s dr/dt dt = int_0^s dr =
                    // r(s) - r(0)
                    real_t r1 = std::sqrt(
                        p.geometry.position.x * p.geometry.position.x +
                        p.geometry.position.y * p.geometry.position.y);
                    real_t r2 = std::sqrt(
                        p_new.geometry.position.x * p_new.geometry.position.x +
                        p_new.geometry.position.y * p_new.geometry.position.y);
                    real_t dr = r2 - r1;
                    // Note, the weight given to this routine already contains
                    // the time step, so we should divide by the time step to
                    // get the average velocity over the time step
                    tally_vector[p.geometry.cell_id] +=
                        dr / (p_new.t - p.t) * weight;
                }
                else {
                    Vec3 vel_tor = transform_cart_to_tor(p.geometry.position,
                                                         p.geometry.velocity);
                    tally_vector[p.geometry.cell_id] += vel_tor.x * weight;
                }
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
                    (0.5 * vel_tor.x * vel_tor.x) * weight;
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
                    (0.5 * vel_tor.y * vel_tor.y) * weight;
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
        case QOI::KDMC_Fluid_Start_Density:
            tally_vector[p.geometry.cell_id] += weight;
            KDMC_avg_fluid_step_duration += p_new.t - p.t;
            KDMC_nb_fluid_steps++;
            break;
        case QOI::KDMC_Fluid_End_Density:
            tally_vector[p_new.geometry.cell_id] += weight;
            KDMC_avg_fluid_step_duration += p_new.t - p.t;
            KDMC_nb_fluid_steps++;
            break;
    }
}

}  // namespace neptune