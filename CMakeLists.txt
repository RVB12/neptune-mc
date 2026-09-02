cmake_minimum_required(VERSION 3.15...3.27)
project(
    neptune
    VERSION 0.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

# ensure -std=c++xx instead of -std=g++xx
set(CMAKE_CXX_EXTENSIONS OFF)

# nicely support folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#
# Get git commit hash for versioning
#
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")
set(NEPTUNE_VERSION "${PROJECT_VERSION}-${GIT_COMMIT_HASH}")

#
# Main library
#
list(APPEND neptune_SOURCES
    src/neptune/boundaryconditions/trim.cpp
    src/neptune/geometry/vec2.cpp
    src/neptune/geometry/vec3.cpp
    src/neptune/geometry/voxelmesh.cpp
    src/neptune/simulator/simulator.cpp
    src/neptune/geometry/trianglemesh.cpp
    src/neptune/geometry/rotatedtrianglemesh.cpp
    src/neptune/geometry/nonuniformmesh.cpp
    src/neptune/geometry/nonuniform1Dmesh.cpp
    src/neptune/geometry/rotatednonuniformmesh.cpp
    src/neptune/physics/update_particle.cpp
    src/neptune/geometry/triangle.cpp
    src/neptune/simulator/simulate_kdmc.cpp
    src/neptune/estimation/estimators.cpp
    src/neptune/simulator/simulate_arcmc.cpp
    src/neptune/estimation/estimator.cpp
    src/neptune/estimation/estimatorgroup.cpp
    src/neptune/sources/pointsource.cpp
    src/neptune/sources/surfacesource.cpp
    src/neptune/sources/volumesource.cpp
    src/neptune/sources/wallreflectionsource.cpp
    src/neptune/sources/recombinationsource.cpp
    src/neptune/sources/amjuelrecombinationsource.cpp
)

add_library(neptune ${neptune_SOURCES})
target_include_directories(neptune PUBLIC src)
target_compile_options(neptune PRIVATE -O3)
target_compile_options(neptune PRIVATE -fPIC)
target_compile_options(neptune PRIVATE -DNDEBUG)
target_compile_definitions(neptune PUBLIC NEPTUNE_VERSION="${NEPTUNE_VERSION}")

include(FetchContent)

# Random123
FetchContent_Declare(
    random123
    GIT_REPOSITORY https://github.com/DEShawResearch/random123.git
    GIT_TAG main
)

FetchContent_MakeAvailable(random123)

add_library(Random123 INTERFACE)
target_include_directories(Random123 INTERFACE
    ${random123_SOURCE_DIR}/include
)

target_link_libraries(neptune PRIVATE Random123)

# JSON

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

target_link_libraries(neptune PUBLIC nlohmann_json::nlohmann_json)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(neptune PUBLIC OpenMP::OpenMP_CXX)
endif()

# PCG

FetchContent_Declare(
    pcg_cpp
    URL https://github.com/imneme/pcg-cpp/archive/428802d1a5634f96bcd0705fab379ff0113bcf13.tar.gz
)
FetchContent_MakeAvailable(pcg_cpp)

add_library(PCGCPP INTERFACE)
target_include_directories(PCGCPP INTERFACE
    ${pcg_cpp_SOURCE_DIR}/include
)

target_link_libraries(neptune PUBLIC PCGCPP)

# Main runner
add_executable(neptune-mc src/neptune/main_hdf5.cpp)
target_link_libraries(neptune-mc PRIVATE neptune)
find_package(HDF5 REQUIRED COMPONENTS CXX)
target_link_libraries(neptune-mc PRIVATE HDF5::HDF5)

#
# Tests
#
find_package(Catch2 3)

if(Catch2_FOUND)
    list(APPEND test_sources
        tests/cpp/test_vec2.cpp
        tests/cpp/test_vec3.cpp
        tests/cpp/test_voxelmesh.cpp
        tests/cpp/test_trianglemesh.cpp
        tests/cpp/test_spatialindex.cpp
        tests/cpp/test_rotatedtrianglemesh.cpp
    )

    add_executable(all_tests ${test_sources})
    target_link_libraries(all_tests PRIVATE Catch2::Catch2WithMain)
    target_link_libraries(all_tests PRIVATE neptune)

    # target_link_libraries(all_tests PRIVATE nlohmann_json::nlohmann_json)
    include(CTest)
    include(Catch)
    catch_discover_tests(all_tests)
endif()