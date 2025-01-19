#ifndef EIRON2_JSON_CONVERSIONS_HPP
#define EIRON2_JSON_CONVERSIONS_HPP

#include <cassert>
#include <nlohmann/json.hpp>

#include "eiron2/common/vec2.hpp"
#include "eiron2/common/vec3.hpp"
#include "eiron2/common/vector.hpp"
#include "eiron2/geometry/voxelmesh.hpp"
#include "eiron2/physics/linesource.hpp"
#include "eiron2/physics/source.hpp"
#include "eiron2/simulator/simulator.hpp"



namespace eiron2 {

void to_json(nlohmann::json& j, const Vec2& v);
void to_json(nlohmann::json& j, const Vec3& v);
void to_json(nlohmann::json& j, const vector<real_t>& v);
void to_json(nlohmann::json& j, const vector<index_t>& v);
void to_json(nlohmann::json& j, const VoxelMesh& m);

void from_json(const nlohmann::json& j, Vec2& v);
void from_json(const nlohmann::json& j, Vec3& v);
void from_json(const nlohmann::json& j, vector<real_t>& v);
void from_json(const nlohmann::json& j, vector<index_t>& v);
void from_json(const nlohmann::json& j, std::unique_ptr<Source>& v);
void from_json(const nlohmann::json& j, VoxelMesh& m);
void from_json(const nlohmann::json& j, std::unique_ptr<Mesh>& v);
void from_json(const nlohmann::json& j, vector<std::unique_ptr<Source>>& v);

}  // namespace eiron2

namespace nlohmann {
template <>
struct adl_serializer<eiron2::PlasmaBackground> {
    static eiron2::PlasmaBackground from_json(const json& j)
    {
        eiron2::index_t n_cells = j.at("plasma_dens").size();
        if(n_cells != (eiron2::index_t) j.at("plasma_temp").size() ||
           n_cells != (eiron2::index_t)  j.at("plasma_vel_x").size() ||
           n_cells != (eiron2::index_t)  j.at("plasma_vel_y").size()){
            throw std::runtime_error("Plasma background arrays must have the same size");
        }

        eiron2::PlasmaBackground m(n_cells);
        j.at("plasma_dens").get_to(m.plasma_dens);
        j.at("plasma_temp").get_to(m.plasma_temp);
        j.at("plasma_vel_x").get_to(m.plasma_vel_x);
        j.at("plasma_vel_y").get_to(m.plasma_vel_y);

        return m;
    }

    static void to_json(json& j, eiron2::PlasmaBackground t)
    {
        j = json{{"plasma_dens", t.plasma_dens},
                 {"plasma_temp", t.plasma_temp},
                 {"plasma_vel_x", t.plasma_vel_x},
                 {"plasma_vel_y", t.plasma_vel_y}};
    }
};

template <>
struct adl_serializer<eiron2::LineSource> {
    static eiron2::LineSource from_json(const json& j)
    {
        eiron2::index_t max_particles = j.at("n_particles");
        eiron2::Vec2 pos_a = j.at("pos_a");
        eiron2::Vec2 pos_b = j.at("pos_b");

        eiron2::LineSource s(max_particles, pos_a, pos_b);

        if (j.find("angle_min") != j.end()) {
            s.angle_min = j.at("angle_min");
        }
        if (j.find("angle_max") != j.end()) {
            s.angle_max = j.at("angle_max");
        }

        return s;
    }

    static void to_json(json& j, eiron2::LineSource t)
    {
        j = json{{"type", "line"},
                 {"pos_a", t.pos_a},
                 {"pos_b", t.pos_b},
                 {"angle_min", t.angle_min},
                 {"angle_max", t.angle_max},
                 {"n_particles", t.max_particles()}};
    }
};

}  // namespace nlohmann

#endif  // EIRON2_JSON_CONVERSIONS_HPP