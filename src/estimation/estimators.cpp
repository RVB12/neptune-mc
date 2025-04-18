#include "eigrene/estimation/estimators.hpp"

namespace eigrene {

    void Estimator::add_track_length(
        Particle& particle,
        Particle& particle_old,
        const Mesh& mesh
    ){
        real_t weight_old = particle_old.weight;
        real_t t_old = particle_old.t;
        index_t cell_id_old = particle_old.geometry.cell_id;
        real_t volume = mesh.get_volume(cell_id_old);

        density[cell_id_old] +=
            (particle.t - t_old) * weight_old / volume * particle.intensity;

        Vec3 m = (particle.t - t_old) * weight_old / volume *
                 particle.intensity * particle.geometry.velocity;
        momentum_x[cell_id_old] += m.x;
        momentum_y[cell_id_old] += m.y;
        momentum_z[cell_id_old] += m.z;

        energy[cell_id_old] +=
            (particle.t - t_old) * weight_old / volume * particle.intensity *
            (0.5 * particle.mass * particle.geometry.velocity.norm2());
    }

    void Estimator::add_final(
        Particle& particle,
        const Mesh& mesh
    ){
        index_t cell_id = particle.geometry.cell_id;
        real_t volume = mesh.get_volume(cell_id);
        density_final[cell_id] +=
            particle.weight / volume * particle.intensity;
    }

    void Estimator::combine(Estimator& other) {
        for (index_t i = 0; i < density.size(); i++) {
            density[i] += other.density[i];
            momentum_x[i] += other.momentum_x[i];
            momentum_y[i] += other.momentum_y[i];
            momentum_z[i] += other.momentum_z[i];
            energy[i] += other.energy[i];
        }
    }

}