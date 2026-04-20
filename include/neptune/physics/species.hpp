#ifndef NEPTUNE_SPECIES_HPP
#define NEPTUNE_SPECIES_HPP
#include <stdexcept>
#include <string>

#include "neptune/common/utils.hpp"

namespace neptune {

const index_t NUM_SPECIES = 3;

enum class Species {
    H = 0,
    D = 1,
    Dplus = 2,
};

inline Species stringtospecies(const std::string& s)
{
    if (s == "H") return Species::H;
    if (s == "D") return Species::D;
    if (s == "D+") return Species::Dplus;
    throw std::invalid_argument("Invalid species string: " + s);
}

inline std::string speciestostring(const Species& species)
{
    switch (species) {
        case Species::H:
            return "H";
        case Species::D:
            return "D";
        case Species::Dplus:
            return "D+";
        default:
            throw std::invalid_argument("Invalid species enum value");
    }
}

// Mass of the species in atomic mass units (u)
const real_t SPECIES_MASS[] = {
    1.00784,  // H
    2.01410,  // D
    2.01410,  // D+
};

}  // namespace neptune

#endif