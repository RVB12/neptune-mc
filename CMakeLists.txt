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
    src/geometry/trianglemesh.cpp
    src/physics/linesource.cpp
    src/physics/source.cpp
    src/physics/update_particle.cpp
    src/common/json_conversions.cpp
)

add_library(eiron2 ${eiron2_SOURCES})
target_include_directories(eiron2 PUBLIC include)
# target_compile_options(eiron2 PRIVATE -O3)
target_compile_options(eiron2 PRIVATE -fPIC)

# JSON
include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

target_link_libraries(eiron2 PUBLIC nlohmann_json::nlohmann_json)


# find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(eiron2 PUBLIC OpenMP::OpenMP_CXX)
# endif()

# Main runner

add_executable(main_runner src/simulator/main.cpp)
target_link_libraries(main_runner PRIVATE eiron2)

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
add_executable(all_tests tests/cpp/test_vec2.cpp tests/cpp/test_vec3.cpp tests/cpp/test_voxelmesh.cpp tests/cpp/test_trianglemesh.cpp tests/cpp/test_json.cpp)
target_link_libraries(all_tests PRIVATE Catch2::Catch2WithMain)
target_link_libraries(all_tests PRIVATE eiron2)
# target_link_libraries(all_tests PRIVATE nlohmann_json::nlohmann_json)

include(CTest)
include(Catch)
catch_discover_tests(all_tests)
# endif()