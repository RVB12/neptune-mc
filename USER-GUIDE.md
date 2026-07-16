# NEPTUNE user guide

This is a user guide for NEPTUNE-MC, we note that documentation for NEPTUNE-FV is not available yet.

## Installation

You should build using CMake. You'll need:

- CMake
- C++ compilers
- Catch2 (This is only needed if you want to run the test suite)
- OpenMP (Optional, but highly recommended to run in parallel)
- HDF5 (Required for the io, easily installed on linux using apt install libhdf5-dev)

To build, run:
- `mkdir build`
- `cmake -B build .`
- `cmake --build build`

This should have made the file: `./build/neptune-mc`.

The executables will need to be on the path, you can do this using:

`export PATH=$PATH:/<path-to-build-folder>`

Note, normally, there would be an install target in CMake that takes care of this for you but I don't know CMake well enough to make that work.

The syntax to use the executable is:

`neptune-mc <input-file>.h5 <output-file>.h5`

If you don't supply an input and output file, the code will use `input.h5` and `output.h5` as defaults.



## The input

You provide input to NEPTUNE in the form of a HDF5 file. HDF5 is a very popular file format (hence, its tools are easily available). The most important advantages over something like JSON/EIRENE input are:
- Availability of tools
- Binary storage (uses much less disk space than simply storing all the numbers in text like in JSON)
- Parallelisability (this is not used right now but could be very useful for domain decomposition in the future)
The disadvantage is that you need to install a separate library

The code expects all input to be present in a single HDF5 file, split into the following groups:
- [the mesh](#the-mesh)
- [the collision processes](#collision-processes)
- [the plasma background](#plasma-background)
- the sources
- the boundary conditions
- the estimators
- other options, like seed, step limit, algorithm, ...

## The Mesh

For the mesh, there are the following options: 
- a structured 3D voxel mesh with uniform cell sizes (`uniform`)
- a structured 3D voxel mesh with nonuniform cell sizes (`non_uniform`)
- an unstructured 2D triangular mesh, where each 3D cell consists of a triangle extruded infinitely in the y-direction (`triangle`)
- an unstructured 2D triangular mesh, where each 3D cell consists of a triangle rotated around the z-axis (`rotated_triangle`)
- a structured 3D voxel mesh with nonuniform cell sizes, where each 3D cell consists of a square rotated around the z-axis (`rotated_non_uniform`)

### Uniform


```JSON
{
    "mesh": {
        "type": "uniform",
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

### Non Uniform

For the structured mesh with nonuniform cell sizes, the input is like so:
```JSON
{
    "mesh": {
        "type": "non_uniform",
        "cb_x": [float], // unit m
        "cb_y": [float], // unit m
        "cb_z": [float], // unit m
    }
}
```
where `cb_x`, `cb_y`, and `cb_z` are arrays of the cell boundaries (in m) in each dimension. If the grid is divided into `nx x n_y x nz` cells, then the dimensions of the arrays are `nx + 1`, `ny + 1`, and `nz + 1` respectively.

### Triangle

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
where `vertices` is a 2D array of vertices (the corners of the triangles). These are 2D points, not 3D! The first column contains the x-dimension of the points and the second column contains the z-dimension of the points. `triangles` is a 2D array of triangle indices with three columns. Each triangle is defined by three vertices.

### Rotated triangle

For the rotated triangle mesh, the input is like so:
```JSON
{
    "mesh": {
        "type": "rotated_triangle",
        "vertices": [[float, float]], // unit m
        "triangles": [[int, int, int]]
    }
}
```
This input is very similar to the input of the triangle mesh, except that first dimension is now the radial direction and the second dimension is the z direction.

### Rotated Non Uniform

For the rotated nonuniform mesh, the input is like so:
```JSON
{
    "mesh": {
        "type": "non_uniform",
        "cb_r": [float], // unit m
        "cb_z": [float], // unit m
    }
}
```
This is similar to the `non_uniform` mesh, except that `cb_r` lists the cell boundaries in the radial dimension and `cb_z` lists the cell boundaries in the z dimension.

## Collision processes

For all collision processes, we assume that the collision rates are expressed w.r.t. time, not distance travelled. The rates are only dependent on the background parameters and do not take the energy of the neutral particle into account. We are currently working on implementing more general collision processes.

For the simple reaction rate formulas, we use rates taken from N. Horsten et al. "Comparison of fluid neutral models for one-
dimensional plasma edge modeling with a finite volume solution of the Boltzmann equation." Physics of Plasmas, 23(1):012510, Jan 2016.

We have the following collision processes:
- `charge_exchange/simple`, charge exchange collisions with rates and post-collisional velocities determined by formulas from Niels' paper.
- `charge_exchange/user_defined`, charge exchange collisions with rates and post-collisional velocities specified per cell.
- `ionisation/simple`, ionisation collisions with rates determined by formulas from Niels' paper.
- `ionisation/user_defined`, ionisation collisions with rates specified per cell.

### Change exchange

A charge exchange collision with simple rate formula is specified as follows

```JSON
{
    "collisions": [
        {
            "type": "charge_exchange",
            "model": "simple",
            "background_species": "H+" or "D+",
            "particle_species": "H" or "D",
        }
    ]
}
```
The formulas for the parameters are:

$R_{cx} = 3.2 \cdot 10^{-15} n_p \sqrt(T_p / 0.026)$, where $R_{cx}$ is the collision rate, $n_p$ is the plasma density and $T_p$ is the plasma temperature.

The post-collisional velocity distribution is a normal distribution with mean the plasma velocity and variance:
$T_p * e / m$, where $e$ converts the temperature from electronvolt to joules and $m$ is the mass of the background particle in kg.

If you wish to specifiy the collision rates yourself, you can do so as follows:

```JSON
{
    "collisions": [
        {
            "type": "charge_exchange",
            "model": "user_defined",
            "background_species": "H+" or "D+",
            "particle_species": "H" or "D",
            "collision_rate": [float], // (in s^{-1})
            "velocity_mean_x": [float], // (in m/s)
            "velocity_mean_y": [float], // (in m/s)
            "velocity_mean_z": [float], // (in m/s)
            "velocity_var": [float],
        }
    ]
}
```
Here, the collision parameters should be defined per cell.

### Ionisation

Ionisation collisions can be specified in a very similar way as charge exchange collisions:


### Change exchange

A charge exchange collision with simple rate formula is specified as follows

```JSON
{
    "collisions": [
        {
            "type": "ionisation",
            "model": "simple",
            "background_species": "H+" or "D+",
            "particle_species": "H" or "D",
        }
    ]
}
```

```JSON
{
    "collisions": [
        {
            "type": "ionisation",
            "model": "user_defined",
            "background_species": "H+" or "D+",
            "particle_species": "H" or "D",
            "collision_rate": [float], // (in s^{-1})
        }
    ]
}
```
The formula for the collision rate is:

$R_{i} = 2.0 \cdot 10^{-13} n_e / (6.0 + T_e / 13.6) \sqrt{T_e / (13.6)} * e^{-13.6 / T_e}$, where $R_{i}$ is the ionisation rate, $n_e$ is the electron density and $T_e$ is the electron temperature.

Notice that the ionisation depends on the electron density and temperature, regardless of the background species that was specified. Don't forget to supply electron parameters if using this model!

## Plasma background

## Sources

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
- `algorithm`: string, either `standard`, `kdmc`, or `arcmc`. By default, the standard algorithm is used.
- `dt`: float, timestep for KDMC. If `kdmc` is selected for the algorithm, this parameter is required.
- `M`: int, parameter for ARC-MC that determines the number of collisions to combine. If `arcmc` is selected, this parameter is required.

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


