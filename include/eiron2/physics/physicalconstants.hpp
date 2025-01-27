#ifndef EIRON2_CONSTANTS_HPP
#define EIRON2_CONSTANTS_HPP

#include <math.h>

#include "eiron2/common/utils.hpp"

namespace eiron2 {
namespace constants {

    // Ion temperature (eV)
    const real_t T_i = 40;
    // Electron temperature (eV)
    const real_t T_e = 40;
    // Charge exchange rate (m^3/s)
    const real_t K_cx = 1.3e-13;
    // Ionization rate (m^3/s)
    const real_t K_i = 2.7e-14;
    // Electron-ion recombination rate (m^3/s)
    const real_t K_r = 1.3e-14;
    // Proton mass (kg)
    const real_t m_p = 1.67e-27;

}  // namespace constants

}  // namespace eiron2
#endif