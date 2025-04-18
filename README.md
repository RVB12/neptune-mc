# EIGRENE

This package simulates neutral particles in a fusion reactor. The name is energie backwards.

### Installation

You should build using CMake. You'll need:

- CMake
- C++ compilers
- Catch2 (This is only needed if you want to run the test suite)

To build, run:
- `mkdir build`
- `cmake -B build .`
- `cmake --build build`

This should have made a file `./build/eigrene-json`. This will read a json file with input parameters and write a json file with the output. Detailed documentation for the json file is pending, but you can have a look at the examples folder in the meantime.

### General design

#### Geometry

I really want to have a strong separation between the geometry and the physics. The idea is that the geometry should only contain geometry things, most importantly, the geometry should be able to answer:
- Is a certain point inside of a cell?
- What cell is a certain point in?
- What is the distance to the nearest boundary?
- What are the neighbours of a cell?
- What is the volume of a cell?
- What is the normal of a boundary?

A downside of this generality is that there can be no optimization for special cases. For example, the boundary conditions use the surface normal for reflections, but this could be simplified for square domains. I believe that trade-off is worth it (we don't care about performance for square domains that much anyway).

Ideally, you can just implement some fancy method on a square domain, and it'll just work when you go to a fancier triangular grid

#### Parallelisation

The main loop is parallelized using OpenMP with private grids, this means that each thread gets its own copy of the density, momentum,... vectors. This reduces communication overhead, but it does of course require a lot of memory for large grids. Should be fine for most use cases though.

Also note that the background parameters like temperature and the definition of the grids are not copied, so that already saves some space over a normal MPI parallelisation.

### Features

The following features are currently supported

- Rectangular grid (structured)
- Triangular mesh grid in x-z direction, extended in y direction with straight lines
- Collision rate and BGK parameters based on background temperature, momentum and density
- time averaged density estimate
- density at a certain time
- KD (only final density)
- Line source
- Point source

The following features are planned

- momentum and energy estimator
- fluid estimator
- Volume source
- Triangular mesh, rotated around z axis to form cells confined by toroidal shaped
- Magnetic sheath BC for the sources
- ...
