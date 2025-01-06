cmake_minimum_required(VERSION 3.15...3.27)
project(
    eiron2
    VERSION 0.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

# ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_EXTENSIONS OFF)

# nicely support folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#
# Main library
#
list(APPEND eiron2_SOURCES
    src/geometry/vec2.cpp
    src/geometry/vec3.cpp
    src/geometry/voxelmesh.cpp
    src/simulator/simulator.cpp
)

add_library(eiron2 ${eiron2_SOURCES})
target_include_directories(eiron2 PUBLIC include)
target_compile_options(eiron2 PRIVATE -O3)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(eiron2 PUBLIC OpenMP::OpenMP_CXX)
endif()

#
# Python bindings
#
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

python_add_library(_core MODULE src/pythonbindings/main.cpp WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers)
target_link_libraries(_core PRIVATE eiron2)

# Note, i just copied this install command from the scikit build example
# this may influence other things, need to check
install(TARGETS _core DESTINATION eiron2)

#
# C++ Examples
#
add_executable(simple_test examples/simple_test.cpp)
target_link_libraries(simple_test PRIVATE eiron2)

# This compile command adds bounds checks for std::vectors,
# but then you also need to compile the library itself with those flags.
# target_compile_options(simple_test PRIVATE -D_GLIBCXX_DEBUG)

#
# Tests
#
find_package(Catch2 3)

# if(Catch2)
add_executable(all_tests tests/cpp/test_vec2.cpp tests/cpp/test_vec3.cpp tests/cpp/test_voxelmesh.cpp)
target_link_libraries(all_tests PRIVATE Catch2::Catch2WithMain)
target_link_libraries(all_tests PRIVATE eiron2)

include(CTest)
include(Catch)
catch_discover_tests(all_tests)
# endif()