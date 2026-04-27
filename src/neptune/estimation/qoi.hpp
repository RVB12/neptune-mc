#ifndef NEPTUNE_QOI_HPP
#define NEPTUNE_QOI_HPP

#include <stdexcept>
#include <string>

namespace neptune {

enum class QOI {
    Density,
    Velocity_x,
    Velocity_y,
    Velocity_z,
    Energy,
    Energy_x,
    Energy_y,
    Energy_z,
};

inline std::string getUnit(QOI qoi)
{
    switch (qoi) {
        case QOI::Density:
            return "m^-3";
        case QOI::Velocity_x:
        case QOI::Velocity_y:
        case QOI::Velocity_z:
            return "(m/s)m^-3";
        case QOI::Energy:
        case QOI::Energy_x:
        case QOI::Energy_y:
        case QOI::Energy_z:
            return "(m/s)^2m^-3";
        default:
            throw std::runtime_error("Unknown QOI");
    }
}

inline QOI stringtoQOI(std::string s)
{
    if (s == "density") return QOI::Density;
    if (s == "velocity_x") return QOI::Velocity_x;
    if (s == "velocity_y") return QOI::Velocity_y;
    if (s == "velocity_z") return QOI::Velocity_z;
    if (s == "energy") return QOI::Energy;
    if (s == "energy_x") return QOI::Energy_x;
    if (s == "energy_y") return QOI::Energy_y;
    if (s == "energy_z") return QOI::Energy_z;
    throw std::runtime_error("Unknown QOI: " + s);
}

inline std::string QOItostring(QOI qoi)
{
    switch (qoi) {
        case QOI::Density:
            return "density";
        case QOI::Velocity_x:
            return "velocity_x";
        case QOI::Velocity_y:
            return "velocity_y";
        case QOI::Velocity_z:
            return "velocity_z";
        case QOI::Energy:
            return "energy";
        case QOI::Energy_x:
            return "energy_x";
        case QOI::Energy_y:
            return "energy_y";
        case QOI::Energy_z:
            return "energy_z";
        default:
            throw std::runtime_error("Unknown QOI");
    }
}
}  // namespace neptune

#endif  // NEPTUNE_QOI_HPP