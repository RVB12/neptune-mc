#include "neptune/common/json_conversions.hpp"

#include <fstream>

namespace neptune {

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


void from_json(const nlohmann::json& j, std::unique_ptr<BoundaryCondition>& m)
{
    std::string type;
    j.at("type").get_to(type);
    if(type == "specular"){
        m = std::make_unique<SpecularBC>();
    } else if(type == "absorbing"){
        m = std::make_unique<AbsorbingBC>();
    } else if(type == "partially_specular"){
        real_t coeff;
        j.at("absorption_coefficient").get_to(coeff);
        m = std::make_unique<PartiallySpecularBC>(coeff);
    } else if(type == "mixed"){
        real_t alpha;
        real_t coeff;
        j.at("alpha").get_to(alpha);
        j.at("absorption_coefficient").get_to(coeff);
        m = std::make_unique<MixedDiffuseReflectiveBC>(alpha, coeff);
    } else {
        throw std::runtime_error("Unknown boundary condition type");
    }
}

void from_json(const nlohmann::json& j, vector<std::unique_ptr<BoundaryCondition>>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

void from_json(const nlohmann::json& j, std::unique_ptr<FluidBoundaryCondition>& m)
{
    std::string type;
    j.at("type").get_to(type);
    if(type == "noflux"){
        m = std::make_unique<NoFluxFluidBoundaryCondition>();
    } else if(type == "absorbing"){
        m = std::make_unique<AbsorbingFluidBoundaryCondition>();
    } else if(type == "robin"){
        real_t alpha = j.at("alpha");
        real_t beta = j.at("beta");
        real_t gamma = j.at("gamma");
        m = std::make_unique<RobinFluidBoundaryCondition>(alpha, beta, gamma);
    } else {
        throw std::runtime_error("Unknown boundary condition type");
    }
}

void from_json(const nlohmann::json& j, vector<std::unique_ptr<FluidBoundaryCondition>>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}



void from_json(const nlohmann::json& j, std::unique_ptr<BGKScatteringProcess>& m)
{

    neptune::index_t n_cells = j.at("rate").size();

    neptune::BGKScatteringProcess bgk(n_cells);
    j.at("rate").get_to(bgk.rate);
    j.at("bgk_stddev").get_to(bgk.bgk_stddev);
    j.at("bgk_mean_x").get_to(bgk.bgk_mean_x);
    j.at("bgk_mean_y").get_to(bgk.bgk_mean_y);
    j.at("bgk_mean_z").get_to(bgk.bgk_mean_z);

    m = std::make_unique<BGKScatteringProcess>(bgk);
}

void from_json(const nlohmann::json& j, vector<std::unique_ptr<BGKScatteringProcess>>& m)
{
    m.resize(j.size());
    for (index_t i = 0; i < m.size(); i++) {
        j.at(i).get_to(m[i]);
    }
}

void from_json(const nlohmann::json& j, std::unique_ptr<IonisationProcess>& m)
{
    neptune::index_t n_cells = j.at("rate").size();

    neptune::IonisationProcess ion(n_cells);
    j.at("rate").get_to(ion.rate);

    m = std::make_unique<IonisationProcess>(ion);
}

void from_json(const nlohmann::json& j, vector<std::unique_ptr<IonisationProcess>>& m)
{
    m.resize(j.size());
    for (index_t i = 0; i < m.size(); i++) {
        j.at(i).get_to(m[i]);
    }
}


void to_json(nlohmann::json& j, const ExecutionStatistics& stats)
{
    j = nlohmann::json{
        {"num_crossings", stats.num_crossings},
        {"num_collisions", stats.num_collisions},
        {"total_particles", stats.total_particles},
        {"grid_building_time", stats.grid_building_time},
        {"trace_time", stats.trace_time},
        {"total_time", stats.total_time},
        {"num_kd_steps", stats.num_kd_steps},
        {"num_kd_to_kinetic_switches", stats.num_kd_to_kinetic_switches},
    };
}

void from_json(const nlohmann::json& j, std::unique_ptr<Source>& s)
{
    std::string type;
    j.at("type").get_to(type);
    if (type == "volume") {
        VolumeSource vs = j;
        s = std::make_unique<VolumeSource>(vs);
    } else if (type == "point") {
        PointSource ps = j;
        s = std::make_unique<PointSource>(ps);
    } else if (type == "surface") {
        SurfaceSource ss = j;
        s = std::make_unique<SurfaceSource>(ss);
    } else {
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
    j.at("Lz").get_to(m.Lz);
    j.at("nx").get_to(m.nx);
    j.at("ny").get_to(m.ny);
    j.at("nz").get_to(m.nz);
    m.dx = m.Lx / m.nx;
    m.dy = m.Ly / m.ny;
    m.dz = m.Lz / m.nz;
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
    } else if (type == "rotated_triangle") {
        RotatedTriangleMesh rtm = j;
        m = std::make_unique<RotatedTriangleMesh>(rtm);
    } else if (type == "nonuniform") {
        NonUniformMesh num = j;
        m = std::make_unique<NonUniformMesh>(num);
    } else if (type == "rotated_nonuniform") {
        RotatedNonUniformMesh num = j;
        m = std::make_unique<RotatedNonUniformMesh>(num);
    }
    else {
        throw std::runtime_error("Unknown mesh type");
    }
}

// void to_json(nlohmann::json& j, const PlasmaBackground& m)
// {
//     j = nlohmann::json{
//         {"plasma_dens", m.plasma_dens},
//         {"plasma_temp", m.plasma_temp},
//         {"plasma_vel_x", m.plasma_vel_x},
//         {"plasma_vel_y", m.plasma_vel_y},
//         {"plasma_vel_z", m.plasma_vel_z},
//     };
// }

// void from_json(const nlohmann::json& j, PlasmaBackground& m)
// {
//     j.at("plasma_dens").get_to(m.plasma_dens);
//     j.at("plasma_temp").get_to(m.plasma_temp);
//     j.at("plasma_vel_x").get_to(m.plasma_vel_x);
//     j.at("plasma_vel_y").get_to(m.plasma_vel_y);
//     j.at("plasma_vel_z").get_to(m.plasma_vel_z);

//     if(j.contains("neutral_dens")) {
//         j.at("neutral_dens").get_to(m.neutral_dens);
//         j.at("neutral_temp").get_to(m.neutral_temp);
//         j.at("neutral_vel_x").get_to(m.neutral_vel_x);
//         j.at("neutral_vel_y").get_to(m.neutral_vel_y);
//         j.at("neutral_vel_z").get_to(m.neutral_vel_z);
//     } else {
//         // set all the values to zero
//         m.neutral_dens.resize(m.plasma_dens.size());
//         m.neutral_temp.resize(m.plasma_dens.size());
//         m.neutral_vel_x.resize(m.plasma_dens.size());
//         m.neutral_vel_y.resize(m.plasma_dens.size());
//         m.neutral_vel_z.resize(m.plasma_dens.size());
//         for(int i = 0; i < m.plasma_dens.size(); i++) {
//             m.neutral_dens[i] = 0;
//             m.neutral_temp[i] = 0;
//             m.neutral_vel_x[i] = 0;
//             m.neutral_vel_y[i] = 0;
//             m.neutral_vel_z[i] = 0;
//         }
//     }

//     if (m.plasma_dens.size() != m.plasma_temp.size() ||
//         m.plasma_dens.size() != m.plasma_vel_x.size() ||
//         m.plasma_dens.size() != m.plasma_vel_y.size() ||
//         m.plasma_dens.size() != m.plasma_vel_z.size()) {
//         throw std::runtime_error(
//             "Plasma background arrays must have the same size");
//     }
// }

void from_json(const nlohmann::json& j, vector<std::unique_ptr<Source>>& v)
{
    v.resize(j.size());
    for (index_t i = 0; i < v.size(); i++) {
        j.at(i).get_to(v[i]);
    }
}

}  // namespace neptune