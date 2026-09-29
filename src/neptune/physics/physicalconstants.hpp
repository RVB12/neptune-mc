#ifndef NEPTUNE_CONSTANTS_HPP
#define NEPTUNE_CONSTANTS_HPP

#include <math.h>

#include "neptune/common/utils.hpp"

namespace neptune {
namespace constants {

    // Proton mass (kg)
    const real_t m_p = 1.67e-27;
    // Elementary charge unit (C)
    const real_t e = 1.6e-19;
    // Joule to eV conversion factor (J/eV)
    const real_t eV_to_Joules = 1.602176634e-19;
    // Atomic mass unit (kg)
    const real_t m_u = 1.660539e-27;

}  // namespace constants

}  // namespace neptune
#endif