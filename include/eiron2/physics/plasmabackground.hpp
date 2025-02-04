#ifndef EIRON2_PLASMA_BACKGROUND
#define EIRON2_PLASMA_BACKGROUND

#include "eiron2/common/utils.hpp"
#include "eiron2/common/vector.hpp"

namespace eiron2 {

class PlasmaBackground {
   public:
    PlasmaBackground() = default;

    PlasmaBackground(index_t n_cells)
        : plasma_dens(n_cells),
          plasma_temp(n_cells),
          plasma_vel_x(n_cells),
          plasma_vel_y(n_cells)
    {}

    // The background plasma density (in m^-3)
    vector<real_t> plasma_dens;
    // The background plasma temperature (in eV)
    vector<real_t> plasma_temp;
    // The background plasma velocity (in m/s)
    vector<real_t> plasma_vel_x;
    // The background plasma velocity (in m/s)
    vector<real_t> plasma_vel_y;
};

}

#endif