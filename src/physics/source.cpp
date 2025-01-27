#include "eiron2/physics/source.hpp"

namespace eiron2 {

index_t Source::num_particles() const { return num_particles_; }

real_t Source::intensity() const { return intensity_; }

}  // namespace eiron2