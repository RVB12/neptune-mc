cmake_minimum_required(VERSION 3.29..3.31)
project(
    eiron2
    VERSION 0.1
    LANGUAGES CXX
)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 20)

    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

list(APPEND eiron2_SOURCES
    src/geometry/vec3.cpp
    src/geometry/voxelmesh.cpp)

# Define the library
add_library(eiron2 ${eiron2_SOURCES})

target_include_directories(eiron2 PUBLIC include)

# The executables
add_executable(simple_test examples/simple_test.cpp)
target_link_libraries(simple_test PRIVATE eiron2)
target_compile_options(simple_test PRIVATE -D_GLIBCXX_DEBUG)


# Tests
find_package(Catch2 3 REQUIRED)
add_subdirectory(tests)

add_executable(all_tests tests/test_vec3.cpp tests/test_voxelmesh.cpp)
target_link_libraries(all_tests PRIVATE Catch2::Catch2WithMain)
target_link_libraries(all_tests PRIVATE eiron2)

include(CTest)
include(Catch)
catch_discover_tests(all_tests)