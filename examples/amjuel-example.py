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
# # 1D test case
#
# In this test case, we implement the third test case from the paper "Hilbert expansion based fluid models for kinetic equations describing
# neutral particles in the plasma edge of a fusion device" V. Maes et al.

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
import json
amjuel_data = json.load(open("../databases/amjuel.json", "r"))
#
# Transform these parameters to a hdf5 file that can be read by the Monte Carlo code
#
with h5py.File(os.path.join("input.h5"), "w") as f:
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
    coll0["model"] = "amjuel"
    coll0["coef_matrix"] = np.array(amjuel_data["lookup"]["H.3"]["3.1.8"][0]["data_blocks"][0]["matrix_t_by_e"]).flatten()
    coll0["output_collision_info"] = 1
    coll1 = collisions.create_group("1")
    coll1["background_species"] = "H+"
    coll1["particle_species"] = "H"
    coll1["type"] = "ionisation"
    coll1["model"] = "amjuel"
    coll1["coef_matrix"] = np.array(amjuel_data["lookup"]["H.4"]["2.1.5"][0]["data_blocks"][0]["matrix_t_by_e"]).flatten()
    coll1["output_collision_info"] = 1
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

# %%
# Run the simulation
start_time = time.time()
output = subprocess.check_output(["neptune-mc"]).decode("utf-8")
end_time = time.time()
print(output)
print(f"Simulation completed in {end_time - start_time:.2f} seconds")

# %%
# Open the output file
with h5py.File("output.h5", "r") as output_file:
    # Plot the QOIs
    H_density = output_file["H"]["track_length"]["density"][:]
    density_unit = output_file["H"]["track_length"]["density"].attrs["unit"]
    density_unit = density_unit.decode("utf-8") if isinstance(density_unit, bytes) else density_unit
    H_velocity_x = output_file["H"]["track_length"]["velocity_x"][:]
    velocity_unit = output_file["H"]["track_length"]["velocity_x"].attrs["unit"]
    velocity_unit = velocity_unit.decode("utf-8") if isinstance(velocity_unit, bytes) else velocity_unit
    H_energy_x = output_file["H"]["track_length"]["energy_x"][:]
    energy_unit = output_file["H"]["track_length"]["energy_x"].attrs["unit"]
    energy_unit = energy_unit.decode("utf-8") if isinstance(energy_unit, bytes) else energy_unit
    plt.figure(figsize=(12, 5))
    plt.subplot(1,3,1)
    plt.plot(xc, H_density, label="H")
    plt.xlabel("x (m)")
    plt.ylabel("Density" + f" ({density_unit})")
    plt.title("Density")
    plt.grid()
    plt.subplot(1,3,2)
    plt.plot(xc, H_velocity_x, label="H")
    plt.xlabel("x (m)")
    plt.ylabel("Velocity" + f" ({velocity_unit})")
    plt.title("Velocity density")
    plt.grid()
    plt.subplot(1,3,3)
    plt.plot(xc, H_energy_x, label="H")
    plt.xlabel("x (m)")
    plt.ylabel("Energy" + f" ({energy_unit})")
    plt.title("Energy density")
    plt.grid()
    plt.tight_layout()

    # Also calculate the temperature in eV to show that the neutral temperature
    # is close to the ion temperature
    H_temperature = (H_energy_x/H_density * 2.0 - (H_velocity_x/H_density )**2) * m / e
    plt.figure()
    plt.plot(xc, H_temperature, label="H")
    plt.plot(xc, ti, label="H+")
    plt.xlabel("x (m)")
    plt.ylabel("Temperature (eV)")
    plt.title("Temperature")
    plt.legend()
    plt.grid()

    plt.show()

# %%
# Compare AMJUEL rates with rates from simple formula

with h5py.File("output.h5", "r") as output_file:
    # Plot the QOIs
    cx_rate_amjuel = output_file["collision_info"]["charge_exchange_0"]["collision_rate"][:]

    K_cx = 3.2 * 1e-15 * np.sqrt(ti / (0.026))
    cx_rate_formula = ni * K_cx

    plt.figure(figsize=(12, 5))
    plt.subplot(1,3,1)
    plt.plot(xc, cx_rate_formula, label="formula rate")
    plt.plot(xc, cx_rate_amjuel, label="amjuel rate")
    plt.xlabel("x (m)")
    plt.ylabel("Collision Rate" + f" ({density_unit})")
    plt.title("R_cx")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    # plt.ylim(1e3, 1e9)

    ion_rate_amjuel = output_file["collision_info"]["ionisation_0"]["collision_rate"][:]
    ne = ni
    te = ti
    K_i = 2.0 * 1e-13 / (6.0 + te / (13.6))\
        * np.sqrt(te / (13.6))\
        * np.exp(-13.6 / (te))
    ion_rate_formula = ne * K_i
    plt.subplot(1,3,2)
    plt.plot(xc, ion_rate_formula, label="formula rate")
    plt.plot(xc, ion_rate_amjuel, label="amjuel rate")
    plt.xlabel("x (m)")
    plt.ylabel("Collision Rate" + f" ({density_unit})")
    plt.title("R_i")
    plt.yscale("log")
    plt.grid()
    plt.legend()

    plt.tight_layout()

    plt.show()
