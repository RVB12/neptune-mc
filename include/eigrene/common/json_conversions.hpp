#ifndef EIGRENE_JSON_CONVERSIONS_HPP
#define EIGRENE_JSON_CONVERSIONS_HPP

#include <cassert>
#include <chrono>
#include <iostream>
#include <nlohmann/json.hpp>

#include "eigrene/common/vector.hpp"
#include "eigrene/geometry/spatialindex.hpp"
#include "eigrene/geometry/trianglemesh.hpp"
#include "eigrene/geometry/vec2.hpp"
#include "eigrene/geometry/vec3.hpp"
#include "eigrene/geometry/voxelmesh.hpp"
#include "eigrene/physics/boundarycondition.hpp"
#include "eigrene/physics/linesource.hpp"
#include "eigrene/physics/pointsource.hpp"
#include "eigrene/physics/source.hpp"
#include "eigrene/simulator/execution_statistics.hpp"
#include "eigrene/simulator/simulator.hpp"

namespace eigrene {

void to_json(nlohmann::json& j, const Vec2& v);
void to_json(nlohmann::json& j, const Vec3& v);
void to_json(nlohmann::json& j, const TriangleIndices& t);
void to_json(nlohmann::json& j, const vector<real_t>& v);
void to_json(nlohmann::json& j, const vector<index_t>& v);
void to_json(nlohmann::json& j, const vector<Vec2>& v);
void to_json(nlohmann::json& j, const vector<TriangleIndices>& v);
void to_json(nlohmann::json& j, const VoxelMesh& m);
void to_json(nlohmann::json& j, const BoundaryCondition& m);
void to_json(nlohmann::json& j, const vector<BoundaryCondition>& v);
void to_json(nlohmann::json& j, const ExecutionStatistics& stats);

void from_json(const nlohmann::json& j, Vec2& v);
void from_json(const nlohmann::json& j, Vec3& v);
void from_json(const nlohmann::json& j, TriangleIndices& t);
void from_json(const nlohmann::json& j, vector<real_t>& v);
void from_json(const nlohmann::json& j, vector<index_t>& v);
void from_json(const nlohmann::json& j, vector<Vec2>& v);
void from_json(const nlohmann::json& j, vector<TriangleIndices>& v);
void from_json(const nlohmann::json& j, std::unique_ptr<Source>& v);
void from_json(const nlohmann::json& j, VoxelMesh& m);
void from_json(const nlohmann::json& j, std::unique_ptr<Mesh>& v);
void from_json(const nlohmann::json& j, vector<std::unique_ptr<Source>>& v);
void from_json(const nlohmann::json& j, BoundaryCondition& v);
void from_json(const nlohmann::json& j, vector<BoundaryCondition>& v);

}  // namespace eigrene

namespace nlohmann {
template <>
struct adl_serializer<eigrene::PlasmaBackground> {
    static eigrene::PlasmaBackground from_json(const json& j)
    {
        eigrene::index_t n_cells = j.at("plasma_dens").size();
        if (n_cells != (eigrene::index_t)j.at("plasma_temp").size() ||
            n_cells != (eigrene::index_t)j.at("plasma_vel_x").size() ||
            n_cells != (eigrene::index_t)j.at("plasma_vel_y").size()) {
            throw std::runtime_error(
                "Plasma background arrays must have the same size");
        }

        eigrene::PlasmaBackground m(n_cells);
        j.at("plasma_dens").get_to(m.plasma_dens);
        j.at("plasma_temp").get_to(m.plasma_temp);
        j.at("plasma_vel_x").get_to(m.plasma_vel_x);
        j.at("plasma_vel_y").get_to(m.plasma_vel_y);

        return m;
    }

    static void to_json(json& j, eigrene::PlasmaBackground t)
    {
        j = json{{"plasma_dens", t.plasma_dens},
                 {"plasma_temp", t.plasma_temp},
                 {"plasma_vel_x", t.plasma_vel_x},
                 {"plasma_vel_y", t.plasma_vel_y}};
    }
};

template <>
struct adl_serializer<eigrene::LineSource> {
    static eigrene::LineSource from_json(const json& j)
    {
        eigrene::index_t num_particles = j.at("n_particles");
        eigrene::real_t intensity = j.at("intensity");
        eigrene::Vec2 pos_a = j.at("pos_a");
        eigrene::Vec2 pos_b = j.at("pos_b");

        eigrene::LineSource s(num_particles, intensity, pos_a, pos_b);

        if (j.find("angle_min") != j.end()) {
            s.angle_min = j.at("angle_min");
        }
        if (j.find("angle_max") != j.end()) {
            s.angle_max = j.at("angle_max");
        }

        return s;
    }

    static void to_json(json& j, eigrene::LineSource t)
    {
        j = json{{"type", "line"},
                 {"pos_a", t.pos_a},
                 {"pos_b", t.pos_b},
                 {"angle_min", t.angle_min},
                 {"angle_max", t.angle_max},
                 {"n_particles", t.num_particles()},
                 {"intensity", t.intensity()}};
    }
};


template <>
struct adl_serializer<eigrene::PointSource> {
    static eigrene::PointSource from_json(const json& j)
    {
        eigrene::index_t num_particles = j.at("n_particles");
        eigrene::real_t intensity = j.at("intensity");
        eigrene::Vec2 pos = j.at("pos");
        eigrene::Vec2 vel_mean = j.at("vel_mean");
        eigrene::real_t vel_stddev = j.at("vel_stddev");

        eigrene::PointSource s(num_particles, intensity, pos, vel_mean, vel_stddev);

        return s;
    }

    static void to_json(json& j, eigrene::PointSource t)
    {
        j = json{{"type", "point"},
                 {"pos", t.pos},
                 {"vel_mean", t.vel_mean},
                 {"vel_stddev", t.vel_stddev},
                 {"n_particles", t.num_particles()},
                 {"intensity", t.intensity()}};
    }
};

template <>
struct adl_serializer<eigrene::TriangleMesh> {
    static eigrene::TriangleMesh from_json(const json& j)
    {
        eigrene::vector<eigrene::Vec2> vertices = j.at("vertices");
        eigrene::vector<eigrene::TriangleIndices> triangles = j.at("triangles");

        return eigrene::TriangleMesh(vertices, triangles);
    }
};

}  // namespace nlohmann

#endif  // EIGRENE_JSON_CONVERSIONS_HPP