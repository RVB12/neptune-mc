#include "neptune/physics/source.hpp"

namespace neptune {

index_t Source::num_particles() const { return num_particles_; }

real_t Source::intensity() const { return intensity_; }

}  // namespace neptune