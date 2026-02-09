#ifndef NEPTUNE_ESTIMATORS_HPP
#define NEPTUNE_ESTIMATORS_HPP

#include "neptune/common/utils.hpp"
#include "neptune/common/vector.hpp"
#include "neptune/finitevolume/fluidboundarycondition.hpp"
#include "neptune/geometry/mesh.hpp"
#include "neptune/physics/particle.hpp"

namespace neptune {

class FluidEstimator {
   public:
    // Constructor
    FluidEstimator(index_t n_cells) : density(n_cells, 0) {}

    FluidEstimator(index_t n_cells, real_t dt_fluid)
        : density(n_cells, 0), dt_fluid(dt_fluid)
    {}

    index_t size() const { return density.size(); }

    /**
     * Add the influence of the particle to the estimates
     * using a track length estimator
     *
     * @param particle     Particle
     *                     The particle
     * @param particle_old Particle
     *                     The particle before the last update
     * @param mesh         Mesh
     *                     The mesh
     */
    void add_fluid_step(const Particle& particle,
                        const Particle& particle_old,
                        const Mesh& mesh);

    /**
     * Combine the estimates of two estimators
     *
     * This is used for parallel simulations
     * where each thread has its own estimator
     * and we want to combine them at the end
     */
    void combine(const FluidEstimator& other);

    /**
     * Do the fluid simulation at the end of the simulation
     */
    void finalize(Mesh& mesh,
                  const vector<std::unique_ptr<FluidBoundaryCondition>>&
                      boundary_conditions,
                  const vector<real_t>& advection_x,
                  const vector<real_t>& advection_y,
                  const vector<real_t>& advection_z,
                  const vector<real_t>& diffusion_rate,
                  const vector<real_t>& diffusion_factor);

    // initial density for fluid solver (and final density afterwards)
    vector<real_t> density;
    // average timestep to be used as final time in the fluid solver
    real_t dt = 0.0;
    // number of particles that were used to estimate the density
    index_t n_particles = 0;
    // Timestep to be used in the fluid solver
    const real_t dt_fluid = 1e-5;
    // Fluid density integral
    vector<real_t> fluid_density_integral;
};

class ArcmcEstimator {
   public:
    // Constructor
    ArcmcEstimator(index_t n_cells)
        : density(n_cells, 0),
          momentum_x(n_cells, 0),
          momentum_y(n_cells, 0),
          momentum_z(n_cells, 0),
          energy(n_cells, 0),
          energy_x(n_cells, 0),
          density_final(n_cells, 0),
          momentum_x_final(n_cells, 0),
          momentum_y_final(n_cells, 0),
          momentum_z_final(n_cells, 0),
          energy_final(n_cells, 0),
          energy_x_final(n_cells, 0)
    {}

    index_t size() const { return density.size(); }

    /**
     * Add the influence of the particle to the estimates
     * using a track length estimator
     *
     * @param particle     Particle
     *                     The particle
     * @param particle_old Particle
     *                     The particle before the last update
     * @param mesh         Mesh
     *                     The mesh
     */
    void add_track_length(const Particle& particle,
                          const Particle& particle_old,
                          const Mesh& mesh,
                          real_t V2_corr,
                          real_t V2_corr_x);
    /**
     * Add the influence of the particle to the estimates
     *
     * This is only meant to be called with particles that
     * satisfy t = t_final, particles that ionise or otherwise
     * dissappear before the final time should not be added
     *
     * @param particle     Particle
     *                     The particle
     * @param mesh         Mesh
     *                     The mesh
     */
    void add_final(const Particle& particle,
                   const Mesh& mesh,
                   real_t V2_corr,
                   real_t V2_corr_x);

    /**
     * Combine the estimates of two estimators
     *
     * This is used for parallel simulations
     * where each thread has its own estimator
     * and we want to combine them at the end
     */
    void combine(const ArcmcEstimator& other);

    // Time averaged density
    vector<real_t> density;
    // Time averaged momentum
    vector<real_t> momentum_x;
    vector<real_t> momentum_y;
    vector<real_t> momentum_z;
    // Time averaged energy
    vector<real_t> energy;
    // Time averaged energy (x only)
    vector<real_t> energy_x;
    // Density at the final time of the simulation
    vector<real_t> density_final;
    // Momentum at the final time of the simulation
    vector<real_t> momentum_x_final;
    vector<real_t> momentum_y_final;
    vector<real_t> momentum_z_final;
    // Energy at the final time of the simulation
    vector<real_t> energy_final;
    vector<real_t> energy_x_final;
};

class Estimator {
   public:
    // Constructor
    Estimator(index_t n_cells)
        : density(n_cells, 0),
          momentum_x(n_cells, 0),
          momentum_y(n_cells, 0),
          momentum_z(n_cells, 0),
          energy(n_cells, 0),
          energy_x(n_cells, 0),
          density_final(n_cells, 0),
          momentum_x_final(n_cells, 0),
          momentum_y_final(n_cells, 0),
          momentum_z_final(n_cells, 0),
          energy_final(n_cells, 0),
          energy_x_final(n_cells, 0)
    {}

    index_t size() const { return density.size(); }

    /**
     * Add the influence of the particle to the estimates
     * using a track length estimator
     *
     * @param particle     Particle
     *                     The particle
     * @param particle_old Particle
     *                     The particle before the last update
     * @param mesh         Mesh
     *                     The mesh
     */
    void add_track_length(const Particle& particle,
                          const Particle& particle_old,
                          const Mesh& mesh);
    /**
     * Add the influence of the particle to the estimates
     *
     * This is only meant to be called with particles that
     * satisfy t = t_final, particles that ionise or otherwise
     * dissappear before the final time should not be added
     *
     * @param particle     Particle
     *                     The particle
     * @param mesh         Mesh
     *                     The mesh
     */
    void add_final(const Particle& particle, const Mesh& mesh);

    /**
     * Combine the estimates of two estimators
     *
     * This is used for parallel simulations
     * where each thread has its own estimator
     * and we want to combine them at the end
     */
    void combine(const Estimator& other);

    void combine(const FluidEstimator& fluid_estimator);

    void combine(const ArcmcEstimator& arcmc_estimator);

    // Time averaged density
    vector<real_t> density;
    // Time averaged momentum
    vector<real_t> momentum_x;
    vector<real_t> momentum_y;
    vector<real_t> momentum_z;
    // Time averaged energy
    vector<real_t> energy;
    // Time averaged energy (x only)
    vector<real_t> energy_x;
    // Density at the final time of the simulation
    vector<real_t> density_final;
    // Momentum at the final time of the simulation
    vector<real_t> momentum_x_final;
    vector<real_t> momentum_y_final;
    vector<real_t> momentum_z_final;
    // Energy at the final time of the simulation
    vector<real_t> energy_final;
    vector<real_t> energy_x_final;
};

#pragma omp declare reduction(EstimatorReduction \
:Estimator : omp_out.combine(omp_in))            \
    initializer(omp_priv = Estimator(omp_orig.size()))

#pragma omp declare reduction(FluidEstimatorReduction \
:FluidEstimator : omp_out.combine(omp_in))            \
    initializer(omp_priv = FluidEstimator(omp_orig.size(), omp_orig.dt_fluid))

#pragma omp declare reduction(ArcmcEstimatorReduction \
:ArcmcEstimator : omp_out.combine(omp_in))            \
    initializer(omp_priv = ArcmcEstimator(omp_orig.size()))

}  // namespace neptune

#endif