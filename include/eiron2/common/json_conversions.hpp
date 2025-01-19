#ifndef EIRON2_JSON_CONVERSIONS_HPP
#define EIRON2_JSON_CONVERSIONS_HPP

#include <cassert>
#include <nlohmann/json.hpp>

#include "eiron2/common/vec2.hpp"
#include "eiron2/common/vec3.hpp"
#include "eiron2/common/vector.hpp"
#include "eiron2/physics/linesource.hpp"
#include "eiron2/physics/source.hpp"
#include "eiron2/geometry/voxelmesh.hpp"
#include "eiron2/simulator/simulator.hpp"

namespace eiron2 {

void to_json(nlohmann::json& j, const Vec2& v);
void to_json(nlohmann::json& j, const Vec3& v);
void to_json(nlohmann::json& j, const vector<real_t>& v);
void to_json(nlohmann::json& j, const vector<index_t>& v);
void to_json(nlohmann::json& j, const LineSource& v);
void to_json(nlohmann::json& j, const VoxelMesh& m);
void to_json(nlohmann::json& j, const PlasmaBackground& m);

void from_json(const nlohmann::json& j, Vec2& v);
void from_json(const nlohmann::json& j, Vec3& v);
void from_json(const nlohmann::json& j, vector<real_t>& v);
void from_json(const nlohmann::json& j, vector<index_t>& v);
void from_json(const nlohmann::json& j, LineSource& v);
void from_json(const nlohmann::json& j, std::unique_ptr<Source>& v);
void from_json(const nlohmann::json& j, VoxelMesh& m);
void from_json(const nlohmann::json& j, std::unique_ptr<Mesh>& v);
void from_json(const nlohmann::json& j, PlasmaBackground& m);
void from_json(const nlohmann::json& j, vector<std::unique_ptr<Source>>& v);

}  // namespace eiron2

#endif  // EIRON2_JSON_CONVERSIONS_HPP