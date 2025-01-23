#ifndef EIRON2_JSON_CONVERSIONS_HPP
#define EIRON2_JSON_CONVERSIONS_HPP

#include <cassert>
#include <nlohmann/json.hpp>
#include <iostream>
#include <chrono>

#include "eiron2/common/vec2.hpp"
#include "eiron2/common/vec3.hpp"
#include "eiron2/common/vector.hpp"
#include "eiron2/geometry/2Dspatialindex.hpp"
#include "eiron2/geometry/trianglemesh.hpp"
#include "eiron2/geometry/voxelmesh.hpp"
#include "eiron2/physics/boundarycondition.hpp"
#include "eiron2/physics/linesource.hpp"
#include "eiron2/physics/source.hpp"
#include "eiron2/simulator/simulator.hpp"

namespace eiron2 {

void to_json(nlohmann::json& j, const Vec2& v);
void to_json(nlohmann::json& j, const Vec3& v);
void to_json(nlohmann::json& j, const Triangle& t);
void to_json(nlohmann::json& j, const vector<real_t>& v);
void to_json(nlohmann::json& j, const vector<index_t>& v);
void to_json(nlohmann::json& j, const vector<Vec2>& v);
void to_json(nlohmann::json& j, const vector<Triangle>& v);
void to_json(nlohmann::json& j, const VoxelMesh& m);
void to_json(nlohmann::json& j, const BoundaryCondition& m);
void to_json(nlohmann::json& j, const vector<BoundaryCondition>& v);

void from_json(const nlohmann::json& j, Vec2& v);
void from_json(const nlohmann::json& j, Vec3& v);
void from_json(const nlohmann::json& j, Triangle& t);
void from_json(const nlohmann::json& j, vector<real_t>& v);
void from_json(const nlohmann::json& j, vector<index_t>& v);
void from_json(const nlohmann::json& j, vector<Vec2>& v);
void from_json(const nlohmann::json& j, vector<Triangle>& v);
void from_json(const nlohmann::json& j, std::unique_ptr<Source>& v);
void from_json(const nlohmann::json& j, VoxelMesh& m);
void from_json(const nlohmann::json& j, std::unique_ptr<Mesh>& v);
void from_json(const nlohmann::json& j, vector<std::unique_ptr<Source>>& v);
void from_json(const nlohmann::json& j, BoundaryCondition& v);
void from_json(const nlohmann::json& j, vector<BoundaryCondition>& v);

}  // namespace eiron2

namespace nlohmann {
template <>
struct adl_serializer<eiron2::PlasmaBackground> {
    static eiron2::PlasmaBackground from_json(const json& j)
    {
        eiron2::index_t n_cells = j.at("plasma_dens").size();
        if (n_cells != (eiron2::index_t)j.at("plasma_temp").size() ||
            n_cells != (eiron2::index_t)j.at("plasma_vel_x").size() ||
            n_cells != (eiron2::index_t)j.at("plasma_vel_y").size()) {
            throw std::runtime_error(
                "Plasma background arrays must have the same size");
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
        eiron2::index_t num_particles = j.at("n_particles");
        eiron2::Vec2 pos_a = j.at("pos_a");
        eiron2::Vec2 pos_b = j.at("pos_b");

        eiron2::LineSource s(num_particles, pos_a, pos_b);

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
        j = json{
            {"type", "line"},           {"pos_a", t.pos_a},
            {"pos_b", t.pos_b},         {"angle_min", t.angle_min},
            {"angle_max", t.angle_max}, {"n_particles", t.num_particles()}};
    }
};

template <>
struct adl_serializer<eiron2::TriangleMesh> {
    static eiron2::TriangleMesh from_json(const json& j)
    {
        eiron2::vector<eiron2::Vec2> vertices = j.at("vertices");
        eiron2::vector<eiron2::Triangle> triangles = j.at("triangles");
        // Hard code nx and ny for now
        eiron2::index_t nx = 50;
        eiron2::index_t ny = 50;
        auto bbox_func = [&vertices, &triangles](eiron2::index_t i) {
            eiron2::Vec2 p1 = vertices[triangles[i].ind1];
            eiron2::Vec2 p2 = vertices[triangles[i].ind2];
            eiron2::Vec2 p3 = vertices[triangles[i].ind3];

            eiron2::real_t min_x = std::min({p1.x, p2.x, p3.x});
            eiron2::real_t max_x = std::max({p1.x, p2.x, p3.x});
            eiron2::real_t min_y = std::min({p1.y, p2.y, p3.y});
            eiron2::real_t max_y = std::max({p1.y, p2.y, p3.y});

            return eiron2::BoundingBox({min_x, min_y}, {max_x, max_y});
        };
        auto time1 = std::chrono::high_resolution_clock::now();
        eiron2::SpatialIndex2D spatialindex(nx, ny, triangles.size(),
                                            bbox_func);
        auto time2 = std::chrono::high_resolution_clock::now();

        eiron2::vector<eiron2::index_t> neighbors;
        eiron2::make_neighborlist(vertices, triangles, spatialindex, neighbors);
        auto time3 = std::chrono::high_resolution_clock::now();

        // eiron2::vector<eiron2::index_t> neighbors2;
        // eiron2::make_neighborlist(triangles, neighbors2);

        // for(int i = 0; i < neighbors.size(); i++){
        //     std::cout << "Neighbor " << i << ": " << neighbors[i] << " " << neighbors2[i] << std::endl;
        // }

        std::cout << "Time to make spatial index: " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << " ms" << std::endl;
        std::cout << "Time to make neighbor list: " << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << " ms" << std::endl;

        return eiron2::TriangleMesh(std::move(vertices), std::move(triangles),
                                    std::move(neighbors));
    }
};

}  // namespace nlohmann

#endif  // EIRON2_JSON_CONVERSIONS_HPP