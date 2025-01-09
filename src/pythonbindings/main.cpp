#include <pybind11/pybind11.h>

#include <eiron2/common/vector.hpp>
#include <eiron2/geometry/mesh.hpp>
#include <eiron2/geometry/voxelmesh.hpp>
#include <eiron2/geometry/trianglemesh.hpp>
#include <eiron2/simulator/simulator.hpp>
#include <string>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) { return i + j; }

namespace py = pybind11;

template <typename T>
void bind_vector(py::module& m, const std::string& typestr)
{
    std::string pyclass_name = std::string("Vector") + typestr;

    py::class_<eiron2::vectorview<T>>(m, pyclass_name.c_str(),
                                      py::buffer_protocol(), py::dynamic_attr())
        // The following lines define a way to convert python classes (like
        // numpy arrays) to the vectorview class
        .def(py::init([](const py::buffer& b) {
            py::buffer_info info = b.request();

            // Check that the python class matches the vector class, i.e. that
            // it has the same type and the correct number of dimensions
            if (info.format != py::format_descriptor<T>::format() ||
                info.ndim != 1) {
                throw std::runtime_error("Incompatible buffer format!");
            }

            // Final step: create the vectorview from the python class and
            // return it.
            auto* v = new eiron2::vectorview<T>(info.shape[0], (T*)info.ptr,
                                                info.strides[0] / sizeof(T));
            return v;
        }));
}

// bind the Cpp mesh class to the python mesh class
void bind_mesh(py::module& m)
{
    py::class_<eiron2::Mesh>(m, "Mesh");

    py::class_<eiron2::VoxelMesh, eiron2::Mesh>(m, "VoxelMesh")
        .def(py::init<double, double, int, int>());

    py::class_<eiron2::TriangleMesh, eiron2::Mesh>(m, "TriangleMesh")
        .def(py::init<eiron2::vectorview<double>, eiron2::vectorview<int>>());
}

// bind the Cpp plasmabackground class to the python plasmabackground class
void bind_background(py::module& m)
{
    py::class_<eiron2::PlasmaBackground>(m, "PlasmaBackground")
        .def(
            py::init<eiron2::vectorview<double>, eiron2::vectorview<double>,
                     eiron2::vectorview<double>, eiron2::vectorview<double>>());
}

void bind_params(py::module& m)
{
    py::class_<eiron2::SimulationParams>(m, "SimulationParams")
        .def(py::init<>())
        .def_readwrite("n_particles", &eiron2::SimulationParams::n_particles)
        .def_readwrite("max_steps", &eiron2::SimulationParams::max_steps)
        .def_readwrite("seed", &eiron2::SimulationParams::seed);
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        Eiron2 - simple C++ monte carlo neutral code
        -----------------------

        .. currentmodule:: eiron2

        .. autosummary::
           :toctree: _generate
    )pbdoc";

    bind_vector<double>(m, "Double");
    bind_vector<int>(m, "Int");
    bind_mesh(m);
    bind_background(m);
    bind_params(m);

    m.def("simulate", &eiron2::simulate, R"pbdoc(
        Simulate particles and tally their paths
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}