#include "eiron2/common/json_conversions.hpp"

#include <fstream>

namespace eiron2 {

void to_json(nlohmann::json& j, const Vec2& v) { j = nlohmann::json{v.x, v.y}; }

void from_json(const nlohmann::json& j, Vec2& v)
{
    j.at(0).get_to(v.x);
    j.at(1).get_to(v.y);
}

void to_json(nlohmann::json& j, const Vec3& v)
{
    j = nlohmann::json{v.x, v.y, v.z};
}

void from_json(const nlohmann::json& j, Vec3& v)
{
    j.at(0).get_to(v.x);
    j.at(1).get_to(v.y);
    j.at(2).get_to(v.z);
}

void to_json(nlohmann::json& j, const TriangleIndices& t)
{
    j = nlohmann::json{t.ind1, t.ind2, t.ind3};
}

void from_json(const nlohmann::json& j, TriangleIndices& t)
{
    j.at(0).get_to(t.ind1);
    j.at(1).get_to(t.ind2);
    j.at(2).get_to(t.ind3);
}

void to_json(nlohmann::json& j, const vector<real_t>& v)
{
    j = nlohmann::json::array();
    for (int i = 0; i < v.size(); i++) {
        j.push_back(v[i]);
    }
}

void from_json(const nlohmann::json& j, vector<real_t>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

void to_json(nlohmann::json& j, const vector<index_t>& v)
{
    j = nlohmann::json::array();
    for (int i = 0; i < v.size(); i++) {
        j.push_back(v[i]);
    }
}

void from_json(const nlohmann::json& j, vector<index_t>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

void to_json(nlohmann::json& j, const vector<Vec2>& v)
{
    j = nlohmann::json::array();
    for (int i = 0; i < v.size(); i++) {
        j.push_back(v[i]);
    }
}

void from_json(const nlohmann::json& j, vector<Vec2>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

void to_json(nlohmann::json& j, const vector<TriangleIndices>& v)
{
    j = nlohmann::json::array();
    for (int i = 0; i < v.size(); i++) {
        j.push_back(v[i]);
    }
}

void from_json(const nlohmann::json& j, vector<TriangleIndices>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

void to_json(nlohmann::json& j, const BoundaryCondition& m)
{
    switch (m) {
        case BoundaryCondition::REFLECTIVE:
            j = "reflective";
            break;
        case BoundaryCondition::ABSORBING:
            j = "absorbing";
            break;
        default:
            throw std::runtime_error("Unknown boundary condition");
    }
}

void from_json(const nlohmann::json& j, BoundaryCondition& m)
{
    std::string s;
    j.get_to(s);
    if (s == "reflective") {
        m = BoundaryCondition::REFLECTIVE;
    }
    else if (s == "absorbing") {
        m = BoundaryCondition::ABSORBING;
    }
    else {
        throw std::runtime_error("Unknown boundary condition");
    }
}

void to_json(nlohmann::json& j, const vector<BoundaryCondition>& v)
{
    j = nlohmann::json::array();
    for (int i = 0; i < v.size(); i++) {
        j.push_back(v[i]);
    }
}

void from_json(const nlohmann::json& j, vector<BoundaryCondition>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

// void to_json(nlohmann::json& j, const LineSource& s)
// {
//     j = nlohmann::json{
//         {"type", "line"},
//         {"pos_a", s.pos_a},
//         {"pos_b", s.pos_b},
//         {"angle_min", s.angle_min},
//         {"angle_max", s.angle_max},
//         {"max_particles", s.max_particles()},
//         {"particle_id_offset", s.particle_id_offset},
//     };
// }

// void from_json(const nlohmann::json& j, LineSource& s)
// {
//     if (j.at("type") != "line") {
//         throw std::runtime_error("Wrong source type");
//     }
//     j.at("pos_a").get_to(s.pos_a);
//     j.at("pos_b").get_to(s.pos_b);
//     j.at("angle_min").get_to(s.angle_min);
//     j.at("angle_max").get_to(s.angle_max);
//     j.at("max_particles").get_to(s.num_particles_);
//     j.at("particle_id_offset").get_to(s.particle_id_offset);
// }

void from_json(const nlohmann::json& j, std::unique_ptr<Source>& s)
{
    std::string type;
    j.at("type").get_to(type);
    if (type == "line") {
        LineSource ls = j;
        s = std::make_unique<LineSource>(ls);
    }
    else {
        throw std::runtime_error("Unknown source type");
    }
}

void to_json(nlohmann::json& j, const VoxelMesh& m)
{
    j = nlohmann::json{
        {"type", "voxel"}, {"Lx", m.Lx}, {"Ly", m.Ly},
        {"nx", m.nx},      {"ny", m.ny},
    };
}

void from_json(const nlohmann::json& j, VoxelMesh& m)
{
    if (j.at("type") != "voxel") {
        throw std::runtime_error("Wrong mesh type");
    }
    j.at("Lx").get_to(m.Lx);
    j.at("Ly").get_to(m.Ly);
    j.at("nx").get_to(m.nx);
    j.at("ny").get_to(m.ny);
    m.dx = m.Lx / m.nx;
    m.dy = m.Ly / m.ny;
}

void from_json(const nlohmann::json& j, std::unique_ptr<Mesh>& m)
{
    std::string type;
    j.at("type").get_to(type);
    if (type == "voxel") {
        VoxelMesh vm = j;
        m = std::make_unique<VoxelMesh>(vm);
    } else if (type == "triangle") {
        TriangleMesh tm = j;
        m = std::make_unique<TriangleMesh>(tm);
    }
    else {
        throw std::runtime_error("Unknown mesh type");
    }
}

void to_json(nlohmann::json& j, const PlasmaBackground& m)
{
    j = nlohmann::json{
        {"plasma_dens", m.plasma_dens},
        {"plasma_temp", m.plasma_temp},
        {"plasma_vel_x", m.plasma_vel_x},
        {"plasma_vel_y", m.plasma_vel_y},
    };
}

void from_json(const nlohmann::json& j, PlasmaBackground& m)
{
    j.at("plasma_dens").get_to(m.plasma_dens);
    j.at("plasma_temp").get_to(m.plasma_temp);
    j.at("plasma_vel_x").get_to(m.plasma_vel_x);
    j.at("plasma_vel_y").get_to(m.plasma_vel_y);

    if (m.plasma_dens.size() != m.plasma_temp.size() ||
        m.plasma_dens.size() != m.plasma_vel_x.size() ||
        m.plasma_dens.size() != m.plasma_vel_y.size()) {
        throw std::runtime_error(
            "Plasma background arrays must have the same size");
    }
}

void from_json(const nlohmann::json& j, vector<std::unique_ptr<Source>>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

}  // namespace eiron2