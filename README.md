# EIGRENE

This package simulates neutral particles in a fusion reactor. The name is energie backwards.

### Installation

#### C++ only

If you only want to use the package from C++, you should build using CMake. You'll need:

- CMake
- C++ compilers
- Catch2 (only if you want to run the tests)

To build, run:
- `mkdir build`
- `cmake -B build .`
- `cmake --build build`

You can now run the simple examples using e.g. `./build/simple_test`.

You can also run `./build/main_runner <input.json> <output.json>`. This will read a json file with input parameters and write a json file with the output. Detailed documentation for the json file is pending, but you can have a look at the examples folder in the meantime.

#### Python bindings

The package can also be used from Python. In essence, this is just a thin wrapper around the C++ code that uses the JSON api. But instead of writing the json files to disk, the string are passed directly, avoiding some file management.

On top of the previous requirements, you'll need:
- Python
- pybind11

To build and install the python bindings, run:
- `pip install .`

You can then use the package from Python as you would expect. For example, you can run the example in the examples folder using `python3 examples/pythontest_2d.py`.

### General design

I really want to have a strong separation between the geometry and the physics. The idea is that the geometry should only contain geometry things, most importantly, the geometry should be able to answer:
- Is a certain point inside of a cell?
- What cell is a certain point in?
- What is the distance to the nearest boundary?
- What are the neighbours of a cell?
- What is the volume of a cell?
- What is the normal of a boundary?

A downside of this generality is that there can be no optimization for special cases. For example, the boundary conditions use the surface normal for reflections, but this could be simplified for square domains. I believe that trade-off is worth it (we don't care about performance for square domains that much anyway).

