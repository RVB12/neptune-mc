#include "eiron2/physics/source.hpp"

namespace eiron2 {

index_t Source::max_particles() const { return max_particles_; }

index_t Source::particles_left() const { return max_particles_ - i_particle; }

}  // namespace eiron2