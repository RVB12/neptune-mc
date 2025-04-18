#ifndef EIGRENE_ESTIMATORS_HPP
#define EIGRENE_ESTIMATORS_HPP

#include "eigrene/common/utils.hpp"
#include "eigrene/common/vector.hpp"
#include "eigrene/geometry/mesh.hpp"
#include "eigrene/physics/particle.hpp"

namespace eigrene {
class Estimator {
   public:
    // Constructor
    Estimator(index_t n_cells)
        : density(n_cells, 0),
          momentum_x(n_cells, 0),
          momentum_y(n_cells, 0),
          momentum_z(n_cells, 0),
          energy(n_cells, 0),
          density_final(n_cells, 0)
    {}

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
    void add_track_length(Particle& particle,
                          Particle& particle_old,
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
    void add_final(Particle& particle, const Mesh& mesh);

    /**
     * Combine the estimates of two estimators
     * 
     * This is used for parallel simulations
     * where each thread has its own estimator
     * and we want to combine them at the end
     */
    void combine(Estimator& other);

    // Time averaged density
    vector<real_t> density;
    // Time averaged momentum
    vector<real_t> momentum_x;
    vector<real_t> momentum_y;
    vector<real_t> momentum_z;
    // Time averaged energy
    vector<real_t> energy;
    // Density at the final time of the simulation
    vector<real_t> density_final;
};
}  // namespace eigrene

#endif