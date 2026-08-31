# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: .venv (3.12.3)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # KDMC convergence example
#
# In this example, we compare KDMC with MC for a homogeneous background with different rates.
#
# You should see that KDMC works very well in the kinetic regime (where it basically revert back to normal MC) and in the diffusive regime.
# The error is large in the intermediate regime.

# %%
import numpy as np
from matplotlib import pyplot as plt
import h5py
import math
import subprocess
import time
import os

# %%

# Problem parameters
nx = 201
ny = 1
nz = 1
Lx = 60.0
Ly = 1.0
Lz = 1.0
n_cells = nx * ny * nz
vol = Lx * Ly * Lz / n_cells

xb = np.linspace(0.0, Lx, nx + 1)
xc = 0.5 * (xb[:-1] + xb[1:])

# manipulate xb so that the first and last cell are massive (to avoid boundary effects)
xb[0] = -100000000.0
xb[-1] = 100000000.0

eps_array = [100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]
eps_array_plots = [100, 1, 0.01]

n_particles = int(1e3)
estimate_variance = True
n_batches = 100
# Note: 1e3*100 = 1e5 particles is not enough to really see the convergence
# 4e7 particles should be better, but that takes a while to run (2 hours on my laptop)
# Note that the variance is only a problem for the diffusive limit, in the kinetic limit, because both method use the same random seed,
# they are actually highly correlated, and the variance is not a problem.
# n_particles = int(4e5)

# %%
#
# Step 2: KDMC simulation
#
# Copy input-mc.h5 to input-kdmc.h5 and modify the parameters for KDMC simulation
for eps in eps_array:
    input_file_mc = f"input-mc-eps{eps}.h5"
    output_file_mc = f"output-mc-eps{eps}.h5"

    R_cx = np.ones(n_cells, dtype=np.float64) / eps
    bgk_var = np.ones(n_cells, dtype=np.float64) / eps
    bgk_mean_x = np.zeros(n_cells, dtype=np.float64)

    with h5py.File(input_file_mc, "w") as f:
        background = f.create_group("background")
        # No need for a background, we will manually define the parameters
        # The code still need a background group, so we create an empty one
        # Mesh
        mesh = f.create_group("mesh")
        mesh["type"] = "non_uniform_1D"
        mesh["cb_x"] = xb
        # Collisions
        collisions = f.create_group("collisions")
        coll0 = collisions.create_group("0")
        coll0["background_species"] = "H+"
        coll0["particle_species"] = "H"
        coll0["type"] = "charge_exchange"
        coll0["model"] = "user_defined"
        coll0["collision_rate"] = R_cx
        coll0["velocity_mean_x"] = bgk_mean_x
        coll0["velocity_mean_y"] = np.zeros_like(bgk_mean_x)
        coll0["velocity_mean_z"] = np.zeros_like(bgk_mean_x)
        coll0["velocity_var"] = bgk_var
        # Sources
        sources = f.create_group("sources")
        sources["n_particles"] = n_particles
        sources["stratification"] = "source_strength"
        # Some volume recombination sources
        source_left = sources.create_group("point_left")
        source_left["type"] = "point"
        source_left["species"] = "H"
        source_left["pos_x"] = 25.0
        source_left["pos_y"] = 0.0
        source_left["pos_z"] = 0.0
        source_left["intensity"] = 1.0
        source_left["vel_mean_x"] = -5.0
        source_left["vel_mean_y"] = 0.0
        source_left["vel_mean_z"] = 0.0
        source_left["vel_std"] = 1.0
        source_right = sources.create_group("point_right")
        source_right["type"] = "point"
        source_right["species"] = "H"
        source_right["pos_x"] = 25.0
        source_right["pos_y"] = 0.0
        source_right["pos_z"] = 0.0
        source_right["intensity"] = 1.0
        source_right["vel_mean_x"] = 15.0
        source_right["vel_mean_y"] = 0.0
        source_right["vel_mean_z"] = 0.0
        source_right["vel_std"] = 1.0
        # Estimators
        estimators = f.create_group("estimators")
        density_estimator = estimators.create_group("density")
        density_estimator["species"] = "H"
        density_estimator["qoi"] = "density"
        density_estimator["estimator_type"] = "final_time"
        density_estimator["estimate_variance"] = 1 if estimate_variance else 0

        # Boundary conditions
        boundary_conditions = f.create_group("boundary_conditions")
        for i in range(6):
            bc = boundary_conditions.create_group(f"bc{i}")
            bc["type"] = "specular"
        # Simulation parameters
        params = f.create_group("params")
        params["max_steps"] = 1000000
        params["n_batches"] = n_batches
        params["seed"] = 1302
        params["t_final"] = 1.0

    # Also create an input file for the KDMC simulation
    input_file_kdmc = f"input-kdmc-eps{eps}.h5"
    output_file_kdmc = f"output-kdmc-eps{eps}.h5"
    with h5py.File(input_file_mc, "r") as f_in:
        with h5py.File(input_file_kdmc, "w") as f_out:
            # Copy all groups and datasets from input-mc.h5 to input-kdmc.h5
            for name in f_in:
                f_in.copy(name, f_out)
            
            # Modify the parameters for KDMC simulation
            f_out["params"]["algorithm"] = "KDMC"
            f_out["params"]["kdmc_dt"] = 1.0

    # Run the MC simulation
    start_time = time.time()
    output = subprocess.check_output(["neptune-mc",input_file_mc ,output_file_mc]).decode("utf-8")
    end_time = time.time()
    mc_time = end_time - start_time

    # Run the KDMC simulation
    start_time = time.time()
    output = subprocess.check_output(["neptune-mc", input_file_kdmc, output_file_kdmc]).decode("utf-8")
    end_time = time.time()
    kdmc_time = end_time - start_time

    print(f"eps={eps}: MC time = {mc_time:.2f} s, KDMC time = {kdmc_time:.2f} s")

    

# %%
# Plot the results

# Plot of the densities for different eps values
for eps in eps_array_plots:
    output_file_mc = f"output-mc-eps{eps}.h5"
    output_file_kdmc = f"output-kdmc-eps{eps}.h5"

    with h5py.File(output_file_mc, "r") as f_mc, h5py.File(output_file_kdmc, "r") as f_kdmc:
        density_mc = f_mc["H"]["final_time"]["density"][:]
        density_kdmc = f_kdmc["H"]["final_time"]["density"][:]

        plt.figure(figsize=(10, 6))
        plt.plot(xc, density_mc, label=f'MC (eps={eps})', linestyle='-')
        plt.plot(xc, density_kdmc, label=f'KDMC (eps={eps})', linestyle="dashed")
        plt.xlabel('x')
        plt.ylabel('Density')
        plt.title(f'Density Comparison for eps={eps}')
        plt.legend()
        plt.grid()

# Also plot the error as a function of eps
errors = []
stddevs = []
for eps in eps_array:
    output_file_mc = f"output-mc-eps{eps}.h5"
    output_file_kdmc = f"output-kdmc-eps{eps}.h5"

    with h5py.File(output_file_mc, "r") as f_mc, h5py.File(output_file_kdmc, "r") as f_kdmc:
        density_mc = f_mc["H"]["final_time"]["density"][:]
        density_kdmc = f_kdmc["H"]["final_time"]["density"][:]

        density_variance_mc = f_mc["H"]["final_time"]["density_variance"][:] / n_batches

        stddevs.append(np.sqrt(np.linalg.norm(density_variance_mc))/np.linalg.norm(density_mc))

        error = np.linalg.norm(density_mc - density_kdmc) / np.linalg.norm(density_mc)
        errors.append(error)

plt.figure(figsize=(10, 6))
plt.loglog(eps_array, errors, marker='o', label='Relative Error')
plt.loglog(eps_array, stddevs, marker='o', label='Standard Deviation')
plt.ylim(1e-4, 1e0)
plt.xlabel('eps')
plt.ylabel('Relative Error')
plt.title('Relative Error between MC and KDMC as a function of eps')
plt.grid()
plt.legend()

# %%
# Clean up the generated files
for eps in eps_array:
    input_file_mc = f"input-mc-eps{eps}.h5"
    output_file_mc = f"output-mc-eps{eps}.h5"
    input_file_kdmc = f"input-kdmc-eps{eps}.h5"
    output_file_kdmc = f"output-kdmc-eps{eps}.h5"

    os.remove(input_file_mc)
    os.remove(output_file_mc)
    os.remove(input_file_kdmc)
    os.remove(output_file_kdmc)
