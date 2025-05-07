#ifndef NEPTUNE_JSON_CONVERSIONS_HPP
#define NEPTUNE_JSON_CONVERSIONS_HPP

#include <cassert>
#include <chrono>
#include <iostream>
#include <nlohmann/json.hpp>

#include "neptune/common/vector.hpp"
#include "neptune/geometry/spatialindex.hpp"
#include "neptune/geometry/trianglemesh.hpp"
#include "neptune/geometry/vec2.hpp"
#include "neptune/geometry/vec3.hpp"
#include "neptune/geometry/voxelmesh.hpp"
#include "neptune/physics/boundarycondition.hpp"
#include "neptune/physics/collision_processes.hpp"
#include "neptune/physics/linesource.hpp"
#include "neptune/physics/pointsource.hpp"
#include "neptune/physics/source.hpp"
#include "neptune/simulator/execution_statistics.hpp"
#include "neptune/simulator/simulator.hpp"

namespace neptune {

void to_json(nlohmann::json& j, const Vec2& v);
void to_json(nlohmann::json& j, const Vec3& v);
void to_json(nlohmann::json& j, const TriangleIndices& t);
void to_json(nlohmann::json& j, const vector<real_t>& v);
void to_json(nlohmann::json& j, const vector<index_t>& v);
void to_json(nlohmann::json& j, const vector<Vec2>& v);
void to_json(nlohmann::json& j, const vector<TriangleIndices>& v);
void to_json(nlohmann::json& j, const VoxelMesh& m);
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
void from_json(const nlohmann::json& j, std::unique_ptr<BoundaryCondition>& v);
void from_json(const nlohmann::json& j,
               vector<std::unique_ptr<BoundaryCondition>>& v);
void from_json(const nlohmann::json& j,
               std::unique_ptr<BGKScatteringProcess>& v);
void from_json(const nlohmann::json& j,
               vector<std::unique_ptr<BGKScatteringProcess>>& m);
void from_json(const nlohmann::json& j, std::unique_ptr<IonisationProcess>& v);
void from_json(const nlohmann::json& j,
               vector<std::unique_ptr<IonisationProcess>>& m);

}  // namespace neptune

namespace nlohmann {
template <>
struct adl_serializer<neptune::PlasmaBackground> {
    static neptune::PlasmaBackground from_json(const json& j)
    {
        neptune::index_t n_cells = j.at("plasma_dens").size();
        if (n_cells != (neptune::index_t)j.at("plasma_temp").size() ||
            n_cells != (neptune::index_t)j.at("plasma_vel_x").size() ||
            n_cells != (neptune::index_t)j.at("plasma_vel_y").size()) {
            throw std::runtime_error(
                "Plasma background arrays must have the same size");
        }

        neptune::PlasmaBackground m(n_cells);
        j.at("plasma_dens").get_to(m.plasma_dens);
        j.at("plasma_temp").get_to(m.plasma_temp);
        j.at("plasma_vel_x").get_to(m.plasma_vel_x);
        j.at("plasma_vel_y").get_to(m.plasma_vel_y);
        j.at("plasma_vel_z").get_to(m.plasma_vel_z);

        if (j.contains("neutral_dens")) {
            j.at("neutral_dens").get_to(m.neutral_dens);
            j.at("neutral_temp").get_to(m.neutral_temp);
            j.at("neutral_vel_x").get_to(m.neutral_vel_x);
            j.at("neutral_vel_y").get_to(m.neutral_vel_y);
            j.at("neutral_vel_z").get_to(m.neutral_vel_z);
        }
        else {
            // set all the values to zero
            m.neutral_dens.resize(m.plasma_dens.size());
            m.neutral_temp.resize(m.plasma_dens.size());
            m.neutral_vel_x.resize(m.plasma_dens.size());
            m.neutral_vel_y.resize(m.plasma_dens.size());
            m.neutral_vel_z.resize(m.plasma_dens.size());
            for (int i = 0; i < m.plasma_dens.size(); i++) {
                m.neutral_dens[i] = 0;
                m.neutral_temp[i] = 0;
                m.neutral_vel_x[i] = 0;
                m.neutral_vel_y[i] = 0;
                m.neutral_vel_z[i] = 0;
            }
        }

        return m;
    }

    static void to_json(json& j, neptune::PlasmaBackground t)
    {
        j = json{{"plasma_dens", t.plasma_dens},
                 {"plasma_temp", t.plasma_temp},
                 {"plasma_vel_x", t.plasma_vel_x},
                 {"plasma_vel_y", t.plasma_vel_y},
                 {"plasma_vel_z", t.plasma_vel_z}};
    }
};

template <>
struct adl_serializer<neptune::LineSource> {
    static neptune::LineSource from_json(const json& j)
    {
        neptune::index_t num_particles = j.at("n_particles");
        neptune::real_t intensity = j.at("intensity");
        neptune::Vec3 pos_a = j.at("pos_a");
        neptune::Vec3 pos_b = j.at("pos_b");
        neptune::Vec3 vel_mean = j.at("vel_mean");
        neptune::real_t vel_stddev = j.at("vel_stddev");

        neptune::LineSource s(num_particles, intensity, pos_a, pos_b, vel_mean,
                              vel_stddev);

        return s;
    }

    static void to_json(json& j, neptune::LineSource t)
    {
        j = json{{"type", "line"},
                 {"pos_a", t.pos_a},
                 {"pos_b", t.pos_b},
                 {"vel_mean", t.vel_mean},
                 {"vel_stddev", t.vel_stddev},
                 {"n_particles", t.num_particles()},
                 {"intensity", t.intensity()}};
    }
};

template <>
struct adl_serializer<neptune::PointSource> {
    static neptune::PointSource from_json(const json& j)
    {
        neptune::index_t num_particles = j.at("n_particles");
        neptune::real_t intensity = j.at("intensity");
        neptune::Vec3 pos = j.at("pos");
        neptune::Vec3 vel_mean = j.at("vel_mean");
        neptune::real_t vel_stddev = j.at("vel_stddev");

        neptune::PointSource s(num_particles, intensity, pos, vel_mean,
                               vel_stddev);

        return s;
    }

    static void to_json(json& j, neptune::PointSource t)
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
struct adl_serializer<neptune::TriangleMesh> {
    static neptune::TriangleMesh from_json(const json& j)
    {
        neptune::vector<neptune::Vec2> vertices = j.at("vertices");
        neptune::vector<neptune::TriangleIndices> triangles = j.at("triangles");

        return neptune::TriangleMesh(vertices, triangles);
    }
};

}  // namespace nlohmann

#endif  // NEPTUNE_JSON_CONVERSIONS_HPP