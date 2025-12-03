# NEPTUNE user guide

This is a user guide for NEPTUNE-MC, we note that documentation for NEPTUNE-FV is not available yet.

### Installation

You should build using CMake. You'll need:

- CMake
- C++ compilers
- Catch2 (This is only needed if you want to run the test suite)
- OpenMP (Optional, but highly recommended to run in parallel)

To build, run:
- `mkdir build`
- `cmake -B build .`
- `cmake --build build`

This should have made the files: `./build/neptune-mc` and `./build/neptune-fv`. These are respectively the Monte-Carlo and finite volume solvers. Both of these programs work by reading a json file as input and the outputting the output in another json file.

The executables will need to be on the path, you can do this using:

`export PATH=$PATH:/<path-to-build-folder>`

Note, normally, there would be an install target in CMake that takes care of this for you but I couldn't get it working.

### The JSON input

You provide input to NEPTUNE in the form of a JSON file. This contains several groups of information:
- the mesh
- the scattering (bgk) processes
- the ionisation processes
- the sources
- the boundary conditions
- other options, like seed, step limit, timestep, ...

#### mesh

For the mesh, there are three options: a structured voxel mesh with uniform cell sizes, a structured voxel mesh with nonuniform cell sizes, or an unstructured triangular mesh.

For the structured mesh with uniform cell sizes, the input is like so:
```JSON
{
    "mesh": {
        "type": "voxel",
        "Lx": float, // unit m
        "Ly": float, // unit m
        "Lz": float, // unit m
        "nx": int,
        "ny": int,
        "nz": int
    }
}
```
where `Lx, Ly, Lz` are the lengths of the domain in each dimension and `nx, ny, nz` are the number of cells in each dimension.

For the structured mesh with nonuniform cell sizes, the input is like so:
```JSON
{
    "mesh": {
        "type": "nonuniform",
        "cb_x": [float], // unit m
        "cb_y": [float], // unit m
        "cb_z": [float], // unit m
    }
}
```
where `cb_x`, `cb_y`, and `cb_z` are arrays of the cell boundaries (in m) in each dimension. If the grid is divided into `nx x n_y x nz` cells, then the dimensions of the arrays are `nx + 1`, `ny + 1`, and `nz + 1` respectively.

For the triangular mesh, the input is like so:
```JSON
{
    "mesh": {
        "type": "triangle",
        "vertices": [[float, float]], // unit m
        "triangles": [[int, int, int]]
    }
}
```
where `vertices` is a list of all the vertices of the mesh (the corners of the triangles). These are 2D points, not 3D! The first dimension of this point is the x-dimension of the full 3D domain and the second dimension of this point is the z-dimension of the full 3D domain. `triangles` is a list of triangles, where each triangle is defined by 3 indices into the vertex list.

#### Collision processes

For all collision processes, we assume that the collision rates are expressed w.r.t. time, not distance travelled. The rates are only dependent on the background parameters and do not take the energy of the neutral particle into account.

For the scattering collision processes, the input is like so:
```JSON
{
    "bgk_processes": [
        {
            "rate": [float], // unit 1/s
            "bgk_stddev": [float], // unit m/s
            "bgk_mean_x": [float], // unit m/s
            "bgk_mean_y": [float], // unit m/s
            "bgk_mean_z": [float], // unit m/s
        }
    ]
}
```
As this is a list, you can define multiple scattering processes this way.

For the ionisation processes, the input is like so:
```JSON
{
    "ionisation_processes": [
        {
            "rate": [float] // unit 1/s
        }
    ]
}
```

All of the vectors (`rate`, `bgk_stddev`, ...) define the rate, bgk_stddev, ... for each cell in the mesh, with the cell index being the index in the rate vector. For example, cell 5 will have rate `rate[5]`.

#### Sources

A point source is defined as follows:

```JSON
{
    "sources": [
        {
            "type": "point",
            "pos": [float, float, float], // unit m
            "vel_mean": [float, float, float], // unit m/s
            "vel_stddev": float, // unit m/s
            "n_particles": int,
            "intensity": float
        }
    ]
}
```

`pos` defines the initial position of the particles. The particles will be initialized with a velocity drawn from the normal distribution with mean `vel_mean` and standard deviation `vel_stddev`. This source uses `n_particles` simulated particles to represent `intensity` physical particles, i.e. each particle from this source represents `intensity/n_particles` particles.

A volume source is defined as follows:

```JSON
{
    "sources": [
        {
            "type": "volume",
            "cell_idx": int,
            "vel_mean": [float, float, float], // unit m/s
            "vel_stddev": float, // unit m/s
            "n_particles": int,
            "intensity": float
        }
    ]
}
```

this is similar to a point source, but instead of on a single point, the particles are generated uniformly within the volume of a specified cell. This is particularly useful for volume recombination.

A surface source is defined as follows:

```JSON
{
    "sources": [
        {
            "type": "surface",
            "cell_idx": int,
            "surface_id": int,
            "vel_magnitude": float, // unit m/s
            "n_particles": int,
            "intensity": float
        }
    ]
}
```

This generates particles on a specified surface of a specified cell. This is particularly useful for recombination at the wall. The velocity is given a specific magnitude (this would typically be the thermal velocity $\sqrt{\frac{8 T}{\pi m}}$) and the direction is sampled from the cosine distribution around the normal.

#### Boundary conditions

Boundary conditions are passed as a list, with the different types shown here:
```JSON
{
    "boundary_conditions": [
        {
            "type": "specular"
        },
        {
            "type": "absorbing"
        },
        {
            "type": "partially_specular",
            "absorption_coefficient": 0.95
        },
    ]
}
```

The order of these boundary conditions shows where they occur in the domain.

For a voxel mesh, the order is: x = 0, y = 0, z = 0, x = Lx, y = Ly, z = Lz. For a triangular mesh, only one boundary condition is supported and it will always select the first one. Allowing multiple boundary conditions in triangular meshes is in the works.

#### Extra parameters

Some extra parameters are:

- `max_steps`: int, the maximum number of steps a single particle can undergo before being removed. This is meant as a protection against infinite loops and we recommend setting this parameter very high to avoid biased simulations.
- `seed`: int, the seed for this simulation. Each particle uses a random number generator with seed `seed + particle_id`. Be careful, this means that simulations with close seeds will be heavily correlated. (Fixing this with proper seed selection for the particles is in the works)
- `t_final`: float, maximum time to simulate the particles. Also enables a special final_time tally
- `dt`: float, timestep for KDMC. If this parameter is set, KDMC will be used instead of MC.

### JSON Output

The JSON output will contain the following arrays:
- `density`: track-length estimate of the particle density (1/m^3)
- `momentum_x`: track-length estimate of the momentum in x-direction (1/(m^2s))
- `momentum_y`: track-length estimate of the momentum in y-direction (1/(m^2s))
- `momentum_z`: track-length estimate of the momentum in z-direction (1/(m^2s))
- `energy`: track-length estimate of the energy
- `density_final`: estimate of the density at t = t_final
- `momentum_x_final`: estimate of the momentum in x-direction at t = t_final
- `momentum_y_final`: estimate of the momentum in y-direction at t = t_final
- `momentum_z_final`: estimate of the momentum in z-direction at t = t_final
- `energy_final`: estimate of the energy at t = t_final

Additionally, a few statistics will be included:
```JSON
{
    "stats": {
        "num_crossings": int, //number of times a particle crossed a cell boundary
        "num_collisions": int, //number of times a particle collided with the plasma background
        "total_particles": int, //total number of particles simulated,
        "grid_building_time": float, //time in seconds spent building the grid, mostly relevant for triangular grid where a spatial index is built up
        "trace_time": float, //time in seconds spent tracing the particles
        "total_time": float, //sum of grid_building_time and trace_time
        "num_kd_steps": int, //number of times a particle did a full kd step
        "num_kd_to_kinetic_switches": int, //number of times a particle needed to switch to kinetic simulation in kd because it hit a boundary

    }
}
```


