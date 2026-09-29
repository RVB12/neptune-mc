# NEPTUNE user guide

This is a user guide for NEPTUNE-MC.

## Installation

You should build using CMake. You'll need:

- CMake
- C++ compilers
- Catch2 (This is only needed if you want to run the test suite)
- OpenMP (Optional, but highly recommended to run in parallel)
- HDF5 (Required for the io, easily installed on linux using apt install libhdf5-dev)

Note, if you are running this on VSC, remember to load the right modules for this, these modules should work:
```
module load foss/2025a
module load matplotlib/3.5.2-foss-2022a
module load GCCcore/14.2.0
module load OpenMPI/5.0.7-GCC-14.2.0
module load CMake/3.31.3-GCCcore-14.2.0
module load Catch2/2.13.10-GCCcore-14.2.0
module load HDF5/1.14.6-gompi-2025a
```
(Thanks Wim for testing this)

The code also depends on a few libraries that will be downloaded automatically by CMake:
- random123
- pcg cpp
- nlohmann json

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
- [the sources](#sources)
- [the boundary conditions](#boundary-conditions)
- [the estimators](#estimators)
- [other options, like seed, step limit, algorithm, ...](#extra-parameters)

The rest of this user guide explains the exact form of this input file. Equally important is the `examples/` folder, which contains a few relevant examples.

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

$R_{cx} = 3.2 \cdot 10^{-15} n_p \sqrt{T_p / 0.026}$, where $R_{cx}$ is the collision rate, $n_p$ is the plasma density and $T_p$ is the plasma temperature.

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

$R_{i} = n_e \frac{2.0 \cdot 10^{-13}}{(6.0 + T_e / 13.6)} \sqrt{\frac{T_e}{13.6}}  e^{-13.6 / T_e}$, where $R_{i}$ is the ionisation rate, $n_e$ is the electron density and $T_e$ is the electron temperature.

Notice that the ionisation depends on the electron density and temperature, regardless of the background species that was specified. Don't forget to supply electron parameters if using this model!

### Outputting collision info

In some cases, you may wish to inspect/use the calculated collision rates and other parameters. The parameters can optionally be added to the output by specifying:
```JSON
{
    "collisions": [
        {
            ...,
            "output_collision_info": 1,
        }
    ]
}
```

## Plasma background

The plasma background is specified as follows:

```JSON
{
    "background": {
        "H+": {
            "density": [float], // unit: m^{-3}
            "velocity_x": [float], // unit: m/s
            "velocity_y": [float], // unit: m/s
            "velocity_z": [float], // unit: m/s
            "temperature": [float], // unit: eV
        },
        "e": {
            ...
        },
        "D+": {
            ...
        },
    },
}
```
The input file must have a `background` group, but you do not need to specify each element, these are optional. The parameters are specified per cell.

## Sources

We support the following particle sources:
- `point`: Particles are generated at a specified point with a velocity drawn from a normal distribution with specified parameters
- `volume`: Particles are generated in a random position within specified cells. The velocities are drawn from a normal distribution with specified parameters
- `recombination`: Particles are generated by a recombination of ions and electrons. Velocities and intensities are calculated using the background parameters using a simple model.
- `wall_recycling`: Particles are generated by ions hitting a boundary and reflecting back as neutrals. Velocities and intensities are calculated using the background parameters using a simple model.

### Number of particles and stratification

Within the sources group, you should also specify how many particles you want to use for the simulation.

```JSON
{
    "sources": {
        "n_particles": int,
        "stratification": "source_strength" or "equal"
    }
}
```

`n_particles` specifies how many particles will be used in total. `stratification` specifies how those particles will distributed among the different sources. With `source_strength`, sources with higher particle intensity will be given more particles. With `equal`, each source will get the same number of particles, regardless of intensity.

### Point source

A point source is defined as follows:

```JSON
{
    "sources": [
        {
            "type": "point",
            "pos_x": float, // unit m
            "pos_y": float, // unit m
            "pos_z": float, // unit m
            "vel_mean_x": float, // unit m/s
            "vel_mean_y": float, // unit m/s
            "vel_mean_z": float, // unit m/s
            "vel_std": float, // unit m/s
            "intensity": float,
            "species": "H" or "D",
        }
    ]
}
```

`pos` defines the initial position of the particles. The particles will be initialized with a velocity drawn from the normal distribution with mean `vel_mean` and standard deviation `vel_std`. This source represents `intensity` physical particles.

### Volume source

A volume source is defined as follows:

```JSON
{
    "sources": [
        {
            "type": "volume",
            "cell_indices": [int],
            "vel_mean_x": [float], // unit m/s
            "vel_mean_y": [float], // unit m/s
            "vel_mean_z": [float], // unit m/s
            "velocity_std": [float], // unit m/s
            "intensities": [float],
            "species": "H" or "D",
        }
    ]
}
```

Notice how the parameters are arrays here. A volume source is often spread over multiple cells. You specify which cells and then specify the parameters per cell (in the same order).

### Recombination source

```JSON
{
    "sources": [
        {
            "type": "recombination",
            "cell_indices": [int],
            "background_species": "H+" or "D+",
            "particle_species": "H" or "D",
        }
    ]
}
```

This is very similar to a volume source, except that the parameters are determined using the plasma background.

The source intensity per cell is:

$0.7 \cdot 10^{-19} \sqrt{\frac{13.6}{T_e}} n_e n_p V$, where $n_e$ is the electron density, $n_p$ is the ion density, $T_e$ is the electron temperature and $V$ is the volume of the cell.

The velocity is drawn from the same Maxwellian as the post-collisional velocity for simple charge exchange.

### Wall Recycling source

```JSON
{
    "sources": [
        {
            "type": "wall_recycling",
            "cell_indices": [int],
            "surface_indices": [int],
            "background_species": "H+" or "D+",
            "particle_species": "H" or "D",
            "sheath_potential": float,
            "franck_condon_temp": float,
            "Rt": float,
            "Rr": float,
        }
    ]
}
```

This is a rather complex source. See Niels' paper for an explanation of the parameters.


## Boundary conditions

For the boundary conditions, you specify a list as follows:
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
        ...
    ]
}
```

Depending on the mesh, different boundary conditions from this list will be applied to different surfaces.

- `uniform`, `non_uniform`:
    - 0: x = 0
    - 1: y = 0
    - 2: z = 0
    - 3: x = Lx
    - 4: y = Ly
    - 5: z = Lz
- `triangle`, `rotated_triangle`:
    - 0: all
    - A more general implementation allowing you to specify which BC is to be applied to which surface is planned

Note that this also implies that each species will undergo the same boundary conditions. We plan to change this in the future.

## Extra parameters

Some extra parameters are:

- `max_steps`: int, the maximum number of steps a single particle can undergo before being removed. This is meant as a protection against infinite loops and we recommend setting this parameter very high to avoid biased simulations. (TODO: add a warning in the output if this is tripped?) (default, 100000)
- `seed`: int, the seed for this simulation. Each particle uses a random number generator with seed `seed_sequence(seed, seed>>32, particle_id, particle_id>>32)`.
- `t_final`: float, maximum time to simulate the particles. (default, -1 -> no limit)
- `n_batches`: number of batches to use (default 1)
- `algorithm`: string, either `MC`, `KDMC`, `KDMC_OLD`, or `ARCMC`. (default, `MC`)
- `kdmc_dt`: float, timestep for KDMC. If `KDMC` or `KDMC_OLD` is selected for the algorithm, this parameter is required.
- `arcmc_M`: int, parameter for ARC-MC that determines the number of collisions to combine. If `ARCMC` is selected, this parameter is required.

## Estimators

Running the code would not be very useful if it didn't generate any output. You can specify which output you want with the estimator group.

```JSON
{
    "estimators": [
        {
            "species": "H" or "D",
            "estimator_type": "track_length" or "final_time",
            "qoi": see below,
            "estimate_variance": 0 or 1
        },
    ]
}
```

### Type

We support two estimator types: `track_length` and `final_time`. The track length estimator keeps track of the time a particle spends in each cell and uses that for estimation. The final time estimator only counts the position of the particle at time `t_final`.

### QOI
Each estimator type can handle several different quantities of interest:
- `density`: The particle density in a cell (m^{-3})
- `velocity_x`: The velocity density in the x direction in a cell, divide by the density to get the average velocity of the particle (m/s m^{-3})
- `velocity_y`: The velocity density in the y direction in a cell, divide by the density to get the average velocity of the particle (m/s m^{-3})
- `velocity_z`: The velocity density in the z direction in a cell, divide by the density to get the average velocity of the particle (m/s m^{-3})
- `energy`: The energy density in a cell, divide by the density to get the average energy of the particle (m^2/s^2 m^{-3})
- `energy_x`: The energy density in a cell, only taking the x direction into account, divide by the density to get the average energy of the particle (m^2/s^2 m^{-3})
- `energy_y`: The energy density in a cell, only taking the y direction into account, divide by the density to get the average energy of the particle (m^2/s^2 m^{-3})
- `energy_z`: The energy density in a cell, only taking the z direction into account, divide by the density to get the average energy of the particle (m^2/s^2 m^{-3})

### Variance
For each of these estimates, you can also ask to get an estimate of the variance. This is useful to estimate how many particles you should be using. This variance is estimated by running the algorithm in multiple batches. The calculated variance is an estimate of the variance of a single batch. Note that the returned estimate of the QOI will be the mean over all the batches, which will have significantly lower variance than a single batch.

An important note is that this can have a meaningful impact on the runtime of the algorithm. If you split 10000 particles into 100 batches of 100 particles, the code can only parallelise over 100 particles at a time. Storing a vector with all the variances also means more memory is used.

### The output

The output of the requested estimators will be written to an output file structured as follows:

```JSON
{
    "D": {
        "track_length": {
            "density": [float],
            "energy": [float],
            "density_variance": [float],
        }
    }
}
```

### KDMC
We also have a special estimator for KDMC. This stores the required input for a fluid estimator: the densities at the start of the diffusive step and the densities at the end of the diffusive step. You can specify these using:

```JSON
{
    "estimators": [
        {
            "species": "H" or "D",
            "estimator_type": "KDMC",
            "qoi": "kdmc_fluid_start_density",
            "estimate_variance": 0 or 1
        },
        {
            "species": "H" or "D",
            "estimator_type": "KDMC",
            "qoi": "kdmc_fluid_end_density",
            "estimate_variance": 0 or 1
        },
    ]
}
```

Getting the final estimate requires some extra post-processing with a fluid solver.