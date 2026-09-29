#ifndef NEPTUNE_RECOMBINATIONSOURCE_HPP
#define NEPTUNE_RECOMBINATIONSOURCE_HPP

#include <cassert>
#include <random>

#include "neptune/common/utils.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/physics/particle.hpp"
#include "neptune/sources/source.hpp"

namespace neptune {

/**
 * @brief A volume source of particles
 *
 * The particles will be generated somewhere in a specified cell.
 *
 * For now, this always models D+ recombining with an electron to form a D atom.
 * Could be extended in the future.
 *
 */
class RecombinationSource : public Source {
   public:
    RecombinationSource(index_t num_particles,
                        index_t cell_idx_,
                        Species background_species_,
                        Species species_,
                        const Mesh& mesh,
                        const Background& background);

    /**
     * @brief Generate a particle
     */
    Particle generate_particle(
        uint64_t particle_id,
        uint64_t seed,
        const Mesh& mesh,
        const Background& background,
        const vector<std::unique_ptr<BoundaryCondition>>& boundary_conditions);

    index_t cell_idx = 0;  // The cell index where the particles are generated
    Species background_species;  // Species of the background the recombines to
                                 // form a neutral
    Species species;             // Species of the generated particle
};

}  // namespace neptune

#endif