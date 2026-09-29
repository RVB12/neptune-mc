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
# # KDMC example
#
# In this test case, we implement the third test case from the paper "Hilbert expansion based fluid models for kinetic equations describing
# neutral particles in the plasma edge of a fusion device" V. Maes et al. This is the same test case as in `1d-example.py`. However, now we also solve it using KDMC.

# %%
import numpy as np
from matplotlib import pyplot as plt
import h5py
import math
import subprocess
import time
import os

# %%

# Physics constants
e = 1.602176634 * 1e-19 # Elementary charge
Da_to_kg = 1.66053906892 * 1e-27 # Dalton (unified atomic mass unit) to kg
m = 1.00784 * Da_to_kg # Mass of hydrogen atom in kg
eV_to_J = 1.602176634 * 1e-19 # Electronvolt to Joules
# Problem parameters
nx = 200
ny = 1
nz = 1
Lx = 1.0
Ly = 1.0
Lz = 1.0
n_cells = nx * ny * nz
vol = Lx * Ly * Lz / n_cells
xb = np.linspace(0, Lx, nx + 1)
xc = 0.5 * (xb[:-1] + xb[1:])
ti = 5.5 + 4.5 * np.cos(2 * math.pi * xc / Lx)
sigma_p_c = np.sqrt(e * ti / m)
gas_puff_dens = 1e21 * np.exp(- 0.5 * (xc - 0.5)**2 / 0.1**2)
ni = 1e21 * np.ones_like(xc)
mean_vel = 0.1 * sigma_p_c

n_particles = 10000
estimate_variance = False
n_batches = 1

# %%
#
# Step 1: normal MC simulation
#
with h5py.File(os.path.join("input-mc.h5"), "w") as f:
    # Background
    background = f.create_group("background")
    Dplus = background.create_group("H+")
    Dplus["density"] = ni
    Dplus["density"].attrs["unit"] = "m^-3"
    Dplus["velocity_x"] = mean_vel
    Dplus["velocity_x"].attrs["unit"] = "m/s"
    Dplus["velocity_y"] = np.zeros_like(xc)
    Dplus["velocity_y"].attrs["unit"] = "m/s"
    Dplus["velocity_z"] = np.zeros_like(xc)
    Dplus["velocity_z"].attrs["unit"] = "m/s"
    Dplus["temperature"] = ti
    Dplus["temperature"].attrs["unit"] = "eV"
    electron = background.create_group("e")
    electron["density"] = ni
    electron["density"].attrs["units"] = "m^-3"
    electron["temperature"] = ti
    electron["temperature"].attrs["units"] = "eV"
    electron["velocity_x"] = np.zeros_like(xc)
    electron["velocity_y"] = np.zeros_like(xc)
    electron["velocity_z"] = np.zeros_like(xc)
    # Mesh
    mesh = f.create_group("mesh")
    mesh["type"] = "uniform"
    mesh["Lx"] = Lx
    mesh["Lx"].attrs["unit"] = "m"
    mesh["Ly"] = Ly
    mesh["Ly"].attrs["unit"] = "m"
    mesh["Lz"] = Lz
    mesh["Lz"].attrs["unit"] = "m"
    mesh["nx"] = nx
    mesh["ny"] = ny
    mesh["nz"] = nz
    # Collisions
    collisions = f.create_group("collisions")
    coll0 = collisions.create_group("0")
    coll0["background_species"] = "H+"
    coll0["particle_species"] = "H"
    coll0["type"] = "charge_exchange"
    coll0["model"] = "simple"
    coll1 = collisions.create_group("1")
    coll1["background_species"] = "H+"
    coll1["particle_species"] = "H"
    coll1["type"] = "ionisation"
    coll1["model"] = "simple"
    # Sources
    sources = f.create_group("sources")
    sources["n_particles"] = n_particles
    sources["stratification"] = "source_strength"
    # Some volume recombination sources
    source = sources.create_group("R_vol")
    source["type"] = "recombination"
    source["background_species"] = "H+"
    source["particle_species"] = "H"
    source["cell_indices"] = np.arange(n_cells, dtype=np.int64)
    source2 = sources.create_group("R_puff")
    source2["type"] = "volume"
    source2["species"] = "H"
    source2["cell_indices"] = np.arange(n_cells, dtype=np.int64)
    source2["intensities"] = gas_puff_dens * vol
    source2["intensities"].attrs["unit"] = "particles/s"
    source2["velocity_x"] = 0.1 * sigma_p_c
    source2["velocity_x"].attrs["unit"] = "m/s"
    source2["velocity_y"] = np.zeros_like(xc)
    source2["velocity_y"].attrs["unit"] = "m/s"
    source2["velocity_z"] = np.zeros_like(xc)
    source2["velocity_z"].attrs["unit"] = "m/s"
    source2["velocity_std"] = sigma_p_c
    source2["velocity_std"].attrs["unit"] = "m/s"
    # Estimators
    estimators = f.create_group("estimators")
    density_estimator = estimators.create_group("density")
    density_estimator["species"] = "H"
    density_estimator["qoi"] = "density"
    density_estimator["estimator_type"] = "track_length"
    density_estimator["estimate_variance"] = 1 if estimate_variance else 0
    velocity_estimator = estimators.create_group("velocity")
    velocity_estimator["species"] = "H"
    velocity_estimator["qoi"] = "velocity_x"
    velocity_estimator["estimator_type"] = "track_length"
    velocity_estimator["estimate_variance"] = 1 if estimate_variance else 0
    energy_estimator = estimators.create_group("energy")
    energy_estimator["species"] = "H"
    energy_estimator["qoi"] = "energy_x"
    energy_estimator["estimator_type"] = "track_length"
    energy_estimator["estimate_variance"] = 1 if estimate_variance else 0
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

# Run the simulation
start_time = time.time()
output = subprocess.check_output(["neptune-mc", "input-mc.h5", "output-mc.h5"]).decode("utf-8")
end_time = time.time()
print(output)
print(f"Simulation completed in {end_time - start_time:.2f} seconds")

# %%
#
# Step 2: KDMC simulation
#
# Copy input-mc.h5 to input-kdmc.h5 and modify the parameters for KDMC simulation
with h5py.File("input-mc.h5", "r") as f_in:
    with h5py.File("input-kdmc.h5", "w") as f_out:
        # Copy all groups and datasets from input-mc.h5 to input-kdmc.h5
        for name in f_in:
            f_in.copy(name, f_out)
        
        # Modify the parameters for KDMC simulation
        f_out["params"]["algorithm"] = "KDMC"
        f_out["params"]["kdmc_dt"] = 1e-6

        # Also add estimators for KDMC fluid start and end density
        estimators = f_out["estimators"]
        kdmc_start_density_estimator = estimators.create_group("kdmc_start_density")
        kdmc_start_density_estimator["species"] = "H"
        kdmc_start_density_estimator["qoi"] = "kdmc_fluid_start_density"
        kdmc_start_density_estimator["estimator_type"] = "KDMC"
        kdmc_start_density_estimator["estimate_variance"] = 1 if estimate_variance else 0
        kdmc_end_density_estimator = estimators.create_group("kdmc_end_density")
        kdmc_end_density_estimator["species"] = "H"
        kdmc_end_density_estimator["qoi"] = "kdmc_fluid_end_density"
        kdmc_end_density_estimator["estimator_type"] = "KDMC"
        kdmc_end_density_estimator["estimate_variance"] = 1 if estimate_variance else 0

        # We will need the collision info to run the fluid solver, so we need to set the output_collision_info flag
        f_out["collisions"]["0"]["output_collision_info"] = 1
        f_out["collisions"]["1"]["output_collision_info"] = 1


# Run the simulation
start_time = time.time()
output = subprocess.check_output(["neptune-mc", "input-kdmc.h5", "output-kdmc.h5"]).decode("utf-8")
end_time = time.time()
print(output)
print(f"Simulation completed in {end_time - start_time:.2f} seconds")

    

# %%
import advection_diffusion
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# Run a fluid solver as post-processing of the KDMC estimator results
with h5py.File("output-kdmc.h5", "r+") as f:

    kdmc_fluid_group = f["H"]["kdmc_fluid"]

    kdmc_start_density = kdmc_fluid_group["kdmc_fluid_start_density"][:]
    kdmc_end_density = kdmc_fluid_group["kdmc_fluid_end_density"][:]

    kdmc_source = kdmc_start_density - kdmc_end_density

    R_cx = f["/collision_info/charge_exchange_0/collision_rate"][:]
    bgk_mean_x = f["/collision_info/charge_exchange_0/post_collision_velocity_mean_x"][:]
    bgk_var = f["/collision_info/charge_exchange_0/post_collision_velocity_variance"][:]

    R_i = f["/collision_info/ionisation_0/collision_rate"][:]

    Dr = 1/(R_cx + R_i)
    Df = bgk_var

    # boundary condition shouldn't really matter so we use something simple
    bcL = [1.0, 0.0, 0.0]
    bcR = [1.0, 0.0, 0.0]

    xb = np.linspace(0, Lx, nx + 1)
    xc = 0.5 * (xb[:-1] + xb[1:])

    def ode(t, u):
        return advection_diffusion.advection_diffusion_dudt(u, bgk_mean_x, Dr, Df, R_i, kdmc_source, xb, bcL, bcR)

    # The nonlinear solver is very sensitive to the initial guess, so we do a few explicit time steps as an initial guess
    start_time = time.time()
    tf = 1.0e-9
    sol = solve_ivp(ode, [0, tf], np.zeros(len(xc)), method='RK45')
    density_start = sol.y[:, -1]
    end_time = time.time()
    print(f"Initial time stepping took {end_time - start_time:.4f} seconds")

    # Now use fsolve to find the steady state solution
    # Note, by not supplying a jacobian, fsolve will use a dense finite difference approximation, which is very expensive for large systems.
    # 1. Fsolve does not support sparse jacobians, so we would need to switch to a custom newton implementation
    # 2. Calculating jacobians is annoying and I didn't want to do it
    start_time = time.time()
    density_fv = fsolve(lambda u: advection_diffusion.advection_diffusion_dudt(u, bgk_mean_x, Dr, Df, R_i, kdmc_source, xb, bcL, bcR), density_start)
    end_time = time.time()
    print(f"fsolve method took {end_time - start_time:.4f} seconds")

    # Calculate momentum and energy density using Vince's formulas
    velocity_x_fv = (density_fv * bgk_mean_x - Dr * np.gradient(density_fv*Df, xc, edge_order=2))
    energy_fv = 0.5 * density_fv * (bgk_var + bgk_mean_x**2) - Dr * np.gradient(bgk_mean_x * density_fv * bgk_var, xc, edge_order=2)

    # Write the results to the HDF5 file
    # (overwrite the existing datasets if they exist, otherwise create them)
    if "kdmc_fluid_density" in kdmc_fluid_group:
        del kdmc_fluid_group["kdmc_fluid_density"]
    kdmc_fluid_group["kdmc_fluid_density"] = density_fv

    if "kdmc_fluid_velocity_x" in kdmc_fluid_group:
        del kdmc_fluid_group["kdmc_fluid_velocity_x"]
    kdmc_fluid_group["kdmc_fluid_velocity_x"] = velocity_x_fv

    if "kdmc_fluid_energy" in kdmc_fluid_group:
        del kdmc_fluid_group["kdmc_fluid_energy"]
    kdmc_fluid_group["kdmc_fluid_energy"] = energy_fv





# %%
#
# Make some plots
#
with h5py.File("output-mc.h5", "r") as output_mc:
    with h5py.File("output-kdmc.h5", "r") as output_kdmc:
        # Plot the QOIs
        H_density = output_mc["H"]["track_length"]["density"][:]
        H_density_kdmc = output_kdmc["H"]["track_length"]["density"][:] + output_kdmc["H"]["kdmc_fluid"]["kdmc_fluid_density"][:]
        density_unit = output_mc["H"]["track_length"]["density"].attrs["unit"]
        density_unit = density_unit.decode("utf-8") if isinstance(density_unit, bytes) else density_unit
        H_velocity_x = output_mc["H"]["track_length"]["velocity_x"][:]
        H_velocity_x_kdmc = output_kdmc["H"]["track_length"]["velocity_x"][:] + output_kdmc["H"]["kdmc_fluid"]["kdmc_fluid_velocity_x"][:]
        velocity_unit = output_mc["H"]["track_length"]["velocity_x"].attrs["unit"]
        velocity_unit = velocity_unit.decode("utf-8") if isinstance(velocity_unit, bytes) else velocity_unit
        H_energy_x = output_mc["H"]["track_length"]["energy_x"][:]
        H_energy_x_kdmc = output_kdmc["H"]["track_length"]["energy_x"][:] + output_kdmc["H"]["kdmc_fluid"]["kdmc_fluid_energy"][:]
        energy_unit = output_mc["H"]["track_length"]["energy_x"].attrs["unit"]
        energy_unit = energy_unit.decode("utf-8") if isinstance(energy_unit, bytes) else energy_unit
        plt.figure(figsize=(12, 5))
        plt.subplot(1,3,1)
        plt.plot(xc, H_density, label="H (MC)")
        plt.plot(xc, H_density_kdmc, label="H (KDMC)", linestyle="dashed")
        plt.xlabel("x (m)")
        plt.ylabel("Density" + f" ({density_unit})")
        plt.title("Density")
        plt.grid()
        plt.legend()
        plt.subplot(1,3,2)
        plt.plot(xc, H_velocity_x, label="H")
        plt.plot(xc, H_velocity_x_kdmc, label="H (KDMC)", linestyle="dashed")
        plt.xlabel("x (m)")
        plt.ylabel("Velocity" + f" ({velocity_unit})")
        plt.title("Velocity density")
        plt.grid()
        plt.subplot(1,3,3)
        plt.plot(xc, H_energy_x, label="H")
        plt.plot(xc, H_energy_x_kdmc, label="H (KDMC)", linestyle="dashed")
        plt.xlabel("x (m)")
        plt.ylabel("Energy" + f" ({energy_unit})")
        plt.title("Energy density")
        plt.grid()
        plt.tight_layout()

    plt.show()
