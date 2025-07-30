# NEPTUNE-MC

NEPTUNE-MC (NEutral Particle Tracing Using New and Efficient Monte Carlo code) is a package that simulates neutral particles in a background of ions.

This package is meant as a step between pure research code and more applied production codes. It specifically tests Kinetic Diffusion Monte Carlo (KDMC), which has previously only been tested in 1D cases, on more complex 3D cases with realistic parameters.

The focus of the software design is on:
1. Simplicity and clean design
2. Speed
3. Extendability

### General design

#### Geometry

I really want to have a strong separation between the geometry and the physics. The idea is that the geometry should only contain geometry things and the physics side should not care about any geometry specifics. The physics should only ask the geometry module the following questions:
- Is a certain point inside of a cell?
- What cell is a certain point in?
- What is the distance to the nearest boundary?
- What are the neighbours of a cell?
- What is the volume of a cell?
- What is the normal of a boundary?
- ...

A downside of this generality is that there can be no optimization for special cases. For example, the boundary conditions use the surface normal for reflections, but this could be simplified for square domains. I believe that trade-off is worth it (we don't care about performance for square domains that much anyway).

Ideally, you can just implement some fancy method on a square domain, and it'll just work when you go to a fancier triangular grid

#### Parallelisation

The main loop is parallelized using OpenMP with private grids, this means that each thread gets its own copy of the density, momentum,... vectors. This reduces communication overhead, but it does of course require a lot of memory for large grids. Anything less than 50.000 cells should be fine.

Also note that the background parameters like temperature and the definition of the grids are not copied, so that already saves some space over a normal MPI parallelisation.

### Features

The following features are currently supported

- Rectangular grid (structured)
- Triangular mesh grid in x-z direction, extended in y direction with straight lines
- Collision rate and BGK parameters based on background temperature, momentum and density
- time averaged estimates for density, momentum, and energy
- estimates at a certain time
- KD (only final density)
- Line source
- Point source
- Fluid estimator

The following features are planned

- Surface source
- Volume source
- Triangular mesh, rotated around z axis to form cells confined by toroidal shaped
- Magnetic sheath BC for the sources
- ...

### Contributors

NEPTUNE is mainly developed by Thijs Steel, with some debugging help from Vince Maes.