#ifndef EIRON2_CONSTANTS_HPP
#define EIRON2_CONSTANTS_HPP

#include <math.h>

#include "eiron2/common/utils.hpp"

namespace eiron2 {
namespace constants {

    // // Ion temperature (eV)
    // const real_t T_i = 40;
    // // Electron temperature (eV)
    // const real_t T_e = 40;
    // // Charge exchange rate (m^3/s)
    // const real_t K_cx = 1.3e-13;
    // // Ionization rate (m^3/s)
    // const real_t K_i = 2.7e-14;
    // // Electron-ion recombination rate (m^3/s)
    // const real_t K_r = 1.3e-14;
    // Proton mass (kg)
    const real_t m_p = 1.67e-27;
    // Elementary charge unit (C)
    const real_t e = 1.6e-19;
    // Joule to eV conversion factor (J/eV)
    const real_t eV_to_Joules = 1.602176634e-19;
    // Atomic mass unit (kg)
    const real_t m_u = 1.660539e-27;

    /**
     * Calculate the charge exchange rate
     * 
     * Formula taken from wouter's thesis
     * 
     * @param T_i Ion temperature (eV)
     */
    inline real_t calculate_K_cx(real_t T_i){
        return 3.2e-15 * std::sqrt(T_i/0.026);
    }

    /**
     * Calculate the ionization rate
     * 
     * Formula taken from wouter's thesis
     * 
     * @param T_e Electron temperature (eV)
     */
    inline real_t calculate_K_i(real_t T_e){
        return 2.0e-13 / (6.0 + T_e/13.6) * std::sqrt(T_e/13.6) * std::exp(-13.6/T_e);
    }

    /**
     * Calculate stddev of the velocity distribution (used for BGK operator)
     * 
     * @param T Temperature (eV)
     * 
     */
    inline real_t calculate_BGK_stddev(real_t T){
        // Assume mass of a single proton
        return std::sqrt(T * eV_to_Joules / m_p);
    }

}  // namespace constants

}  // namespace eiron2
#endif