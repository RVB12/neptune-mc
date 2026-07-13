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

# %%
#
# STEP 1: Load and process the data
#

import numpy as np
import math
import json
import os
import time
import matplotlib.pyplot as plt
import h5py
import subprocess

if(not os.path.exists("fusion_bg_for_MC_uniform_400.npz")):
    raise FileNotFoundError("fusion_bg_for_MC_uniform_400.npz not found. It is not included in the repository. Please request it by email.")
fusion_bg_for_MC_uniform_400 = np.load("fusion_bg_for_MC_uniform_400.npz")

# Number of particles to use in the simulation
num_particles = int(1e5)

# Physics constants
e = 1.602176634 * 1e-19 # Elementary charge
Da_to_kg = 1.66053906892 * 1e-27 # Dalton (unified atomic mass unit) to kg
m = 2.0141017778 * Da_to_kg # Mass of deuterium
eV = 1.602176634 * 1e-19 # Electronvolt

Lx = fusion_bg_for_MC_uniform_400["grid"][-1]
Ly = 1.0
Lz = 1.0
nx = fusion_bg_for_MC_uniform_400["grid"].size - 1
ny = 1
nz = 1

ncells = nx * ny * nz

xb = fusion_bg_for_MC_uniform_400["grid"]
xc = fusion_bg_for_MC_uniform_400["grid_mid"]

te = fusion_bg_for_MC_uniform_400["te"] 
ti = fusion_bg_for_MC_uniform_400["ti"] 

# ti = np.zeros_like(xc) + 20 * eV
# te = ti

ni = fusion_bg_for_MC_uniform_400["ni"]
up = fusion_bg_for_MC_uniform_400["up"]
pit = fusion_bg_for_MC_uniform_400["pit"]
# Assume electron and ion densities are balanced (Zi = 1)
ne = ni

# Recombination (expressed as a density: multiply by volume to get a particle rate!!!)
K_r = 0.7 * 1e-19 * np.sqrt(13.6 / te)
volume_recombination = ni * ne * K_r 

# Charge exchange
K_cx = 3.2 * 1e-15 * np.sqrt(ti / (0.026))
R_cx = ni * K_cx
bgk_stddev = np.sqrt(ti / m )
bgk_mean_x = up * pit

# Ionization
K_i = 2.0 * 1e-13 / (6.0 + te / (13.6))\
    * np.sqrt(te / (13.6))\
    * np.exp(-13.6 / (te))
R_i = ne * K_i

# Other parameters
R_n_t = 0.9 # Neutral reflection coefficient at the target.
R_i_t = 0.9 # Recycling rate of ions to neutrals at the target.
R_r = 0.8 # Fraction of recycled particles that are immediatley back-scattered.

alpha = 0.5
# # Electric potential difference due to the electric field of the sheath.
delta_pot_sh = 3.1


# Number of particles to use in the simulation
estimate_variance = False
n_batches = 1

# Plot these parameters to check
plt.figure(figsize=(12, 8))
plt.subplot(3, 2, 1)
plt.plot(xc, te / eV, label="Te")
plt.plot(xc, ti / eV, label="Ti")
plt.xlabel("x (m)")
plt.ylabel("Temperature (eV)")
plt.title("Temperature")
plt.legend()

plt.subplot(3, 2, 2)
plt.plot(xc, ni, label="ni")
plt.xlabel("x (m)")
plt.ylabel("Density (m^-3)")
plt.title("Ion Density")
plt.legend()

plt.subplot(3, 2, 3)
plt.plot(xc, volume_recombination, label="S")
plt.xlabel("x (m)")
plt.ylabel("Rate (m^-3 s^-1)")
plt.title("Volume Recombination source")
plt.yscale("log")
plt.legend()

plt.subplot(3, 2, 4)
plt.plot(xc, R_cx, label="R_cx")
plt.plot(xc, R_i, label="R_i")
plt.xlabel("x (m)")
plt.ylabel("Rate (s^-1)")
plt.title("Charge Exchange and Ionization Rates")
plt.yscale("log")
plt.legend()

plt.subplot(3, 2, 5)
plt.plot(xc, bgk_mean_x)
plt.xlabel("x (m)")
plt.ylabel("v (m/s)")
plt.title("Mean velocity")

plt.subplot(3, 2, 6)
plt.plot(xc, bgk_stddev)
plt.xlabel("x (m)")
plt.ylabel("v (m/s)")
plt.title("Velocity stddev")

plt.tight_layout()

# %%
with h5py.File(os.path.join("input.h5"), "w") as f:
    # Background
    background = f.create_group("background")
    Dplus = background.create_group("D+")
    Dplus["density"] = ni
    Dplus["density"].attrs["unit"] = "m^-3"
    Dplus["velocity_x"] = bgk_mean_x
    Dplus["velocity_x"].attrs["unit"] = "m/s"
    Dplus["velocity_y"] = np.zeros_like(xc)
    Dplus["velocity_y"].attrs["unit"] = "m/s"
    Dplus["velocity_z"] = np.zeros_like(xc)
    Dplus["velocity_z"].attrs["unit"] = "m/s"
    Dplus["temperature"] = ti
    Dplus["temperature"].attrs["unit"] = "eV"
    # 
    # Mesh
    # 
    mesh = f.create_group("mesh")
    mesh["type"] = "non_uniform"
    mesh["cb_x"] = xb
    mesh["cb_y"] = [-Ly/2, Ly/2]
    mesh["cb_z"] = [-Lz/2, Lz/2]
    # 
    # Collisions
    # 
    collisions = f.create_group("collisions")
    coll0 = collisions.create_group("0")
    coll0["background_species"] = "D+"
    coll0["particle_species"] = "D"
    coll0["type"] = "charge_exchange"
    coll0["model"] = "simple"
    coll1 = collisions.create_group("1")
    coll1["background_species"] = "D+"
    coll1["particle_species"] = "D"
    coll1["type"] = "ionisation"
    coll1["model"] = "simple"
    # 
    # Sources
    # 
    sources = f.create_group("sources")
    sources["n_particles"] = num_particles
    sources["stratification"] = "equal"
    
    recombination_source = sources.create_group("R_vol")
    recombination_source["type"] = "recombination"
    recombination_source["background_species"] = "D+"
    recombination_source["particle_species"] = "D"
    recombination_source["cell_indices"] = np.arange(ncells, dtype=np.int64)
    wall_source = sources.create_group("R_wall")
    wall_source["type"] = "wall_recycling"
    wall_source["cell_indices"] = [nx-1]
    wall_source["surface_indices"] = [3]  # The right wall is surface 3
    wall_source["sheath_potential"] = delta_pot_sh
    wall_source["franck_condon_temp"] = 2.0
    wall_source["Rt"] = 0.9
    wall_source["Rr"] = 0.8
    wall_source["background_species"] = "D+"
    wall_source["particle_species"] = "D"

    # 
    # Estimators
    # 
    estimators = f.create_group("estimators")
    density_estimator = estimators.create_group("density")
    density_estimator["species"] = "D"
    density_estimator["qoi"] = "density"
    density_estimator["estimator_type"] = "track_length"
    density_estimator["estimate_variance"] = 1 if estimate_variance else 0
    velocity_estimator = estimators.create_group("velocity")
    velocity_estimator["species"] = "D"
    velocity_estimator["qoi"] = "velocity_x"
    velocity_estimator["estimator_type"] = "track_length"
    velocity_estimator["estimate_variance"] = 1 if estimate_variance else 0
    energy_estimator = estimators.create_group("energy")
    energy_estimator["species"] = "D"
    energy_estimator["qoi"] = "energy_x"
    energy_estimator["estimator_type"] = "track_length"
    energy_estimator["estimate_variance"] = 1 if estimate_variance else 0
    # Boundary conditions
    boundary_conditions = f.create_group("boundary_conditions")
    for i in range(6):
        if i == 0:
            bc = boundary_conditions.create_group(f"bc{i}")
            bc["type"] = "absorbing"
        elif i == 3:
            bc = boundary_conditions.create_group(f"bc{i}")
            bc["type"] = "partially_specular"
            bc["absorption_coefficient"] = 1.0-R_n_t
        else:
            bc = boundary_conditions.create_group(f"bc{i}")
            bc["type"] = "specular"
    # Simulation parameters
    params = f.create_group("params")
    params["max_steps"] = 10000
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
    D_density = output_file["D"]["track_length"]["density"][:]
    density_unit = output_file["D"]["track_length"]["density"].attrs["unit"]
    density_unit = density_unit.decode("utf-8") if isinstance(density_unit, bytes) else density_unit
    D_velocity_x = output_file["D"]["track_length"]["velocity_x"][:]
    velocity_unit = output_file["D"]["track_length"]["velocity_x"].attrs["unit"]
    velocity_unit = velocity_unit.decode("utf-8") if isinstance(velocity_unit, bytes) else velocity_unit
    D_energy_x = output_file["D"]["track_length"]["energy_x"][:]
    energy_unit = output_file["D"]["track_length"]["energy_x"].attrs["unit"]
    energy_unit = energy_unit.decode("utf-8") if isinstance(energy_unit, bytes) else energy_unit
    plt.figure(figsize=(12, 5))
    plt.subplot(1,3,1)
    plt.plot(xc, D_density, label="D")
    plt.xlabel("x (m)")
    plt.ylabel("Density" + f" ({density_unit})")
    plt.title("Density")
    plt.yscale("log")
    plt.grid()

    plt.subplot(1,3,2)
    plt.plot(xc, D_velocity_x/D_density, label="D")
    plt.xlabel("x (m)")
    plt.ylabel("Velocity" + f" ({velocity_unit})")
    plt.title("Mean Velocity")
    plt.ylim(-1000,1000)
    plt.grid()
    
    D_temperature = (D_energy_x/D_density * 2.0 - (D_velocity_x/D_density )**2) * m / e
    plt.subplot(1,3,3)
    plt.plot(xc, D_temperature, label="D")
    plt.plot(xc, ti, label="D+")
    plt.xlabel("x (m)")
    plt.ylabel("eV")
    plt.title("Temperature")
    plt.grid()
    plt.tight_layout()

    plt.show()
