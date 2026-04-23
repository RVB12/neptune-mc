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
#
# We note that in that test case, hydrogen atoms are simulated, whereas we simulate deuterium in this test case.

# %%
import numpy as np
from matplotlib import pyplot as plt
import h5py
import math
import subprocess
import time
import os

# %%
f = h5py.File(os.path.join("input.h5"), "w")
# Physics constants
e = 1.602176634 * 1e-19 # Elementary charge
Da_to_kg = 1.66053906892 * 1e-27 # Dalton (unified atomic mass unit) to kg
m = 1.007825031898 * Da_to_kg # Mass of Hydrogen-1
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

n_particles = 20000
estimate_variance = False
n_batches = 1

# %%
#
# Transform these parameters to a hdf5 file that can be read by the Monte Carlo code
#

# Background
background = f.create_group("background")
Dplus = background.create_group("D+")
Dplus["density"] = ni
Dplus["density"].attrs["units"] = "m^-3"
Dplus["velocity_x"] = mean_vel
Dplus["velocity_x"].attrs["units"] = "m/s"
Dplus["velocity_y"] = np.zeros_like(xc)
Dplus["velocity_y"].attrs["units"] = "m/s"
Dplus["velocity_z"] = np.zeros_like(xc)
Dplus["velocity_z"].attrs["units"] = "m/s"
Dplus["temperature"] = ti
Dplus["temperature"].attrs["units"] = "eV"
# Mesh
mesh = f.create_group("mesh")
mesh["type"] = "uniform"
mesh["Lx"] = Lx
mesh["Lx"].attrs["units"] = "m"
mesh["Ly"] = Ly
mesh["Ly"].attrs["units"] = "m"
mesh["Lz"] = Lz
mesh["Lz"].attrs["units"] = "m"
mesh["nx"] = nx
mesh["ny"] = ny
mesh["nz"] = nz
# Collisions
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
# Sources
sources = f.create_group("sources")
sources["n_particles"] = n_particles
sources["stratification"] = "source_strength"
# Some volume recombination sources
source = sources.create_group("R_vol")
source["type"] = "recombination"
source["background_species"] = "D+"
source["particle_species"] = "D"
source["cell_indices"] = np.arange(n_cells, dtype=np.int64)
source2 = sources.create_group("R_puff")
source2["type"] = "volume"
source2["species"] = "D"
source2["cell_indices"] = np.arange(n_cells, dtype=np.int64)
source2["intensities"] = gas_puff_dens * vol
source2["intensities"].attrs["units"] = "particles/s"
source2["velocity_x"] = 0.1 * sigma_p_c
source2["velocity_x"].attrs["units"] = "m/s"
source2["velocity_y"] = np.zeros_like(xc)
source2["velocity_y"].attrs["units"] = "m/s"
source2["velocity_z"] = np.zeros_like(xc)
source2["velocity_z"].attrs["units"] = "m/s"
source2["velocity_std"] = sigma_p_c
source2["velocity_std"].attrs["units"] = "m/s"
# Estimators
estimators = f.create_group("estimators")
density_estimator = estimators.create_group("density")
density_estimator["species"] = "D"
density_estimator["qoi"] = "density"
density_estimator["estimator_type"] = "track_length"
density_estimator["estimate_variance"] = 1 if estimate_variance else 0
momentum_estimator = estimators.create_group("momentum")
momentum_estimator["species"] = "D"
momentum_estimator["qoi"] = "momentum_x"
momentum_estimator["estimator_type"] = "track_length"
momentum_estimator["estimate_variance"] = 1 if estimate_variance else 0
energy_estimator = estimators.create_group("energy")
energy_estimator["species"] = "D"
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
# Close the file (important, otherwise neptune-mc won't be able to read it)
f.close()

# %%
# Run the simulation
start_time = time.time()
output = subprocess.check_output(["neptune-mc"]).decode("utf-8")
end_time = time.time()
print(output)
print(f"Simulation completed in {end_time - start_time:.2f} seconds")

# %%
# Open the output file
output_file = h5py.File("output.h5", "r")
# Plot the QOIs
D_density = output_file["D"]["track_length"]["density"][:]
D_momentum_x = output_file["D"]["track_length"]["momentum_x"][:]
D_energy_x = output_file["D"]["track_length"]["energy_x"][:]
plt.figure(figsize=(12, 5))
plt.subplot(1,3,1)
plt.plot(xc, D_density, label="D")
plt.xlabel("x (m)")
plt.ylabel("Density (particles/m^3)")
plt.title("Density")
plt.grid()
plt.subplot(1,3,2)
plt.plot(xc, D_momentum_x, label="D")
plt.xlabel("x (m)")
plt.ylabel("Momentum (m/s)/m^3")
plt.title("Momentum density")
plt.grid()
plt.subplot(1,3,3)
plt.plot(xc, D_energy_x, label="D")
plt.xlabel("x (m)")
plt.ylabel("Energy (m^2/s^2)/m^3")
plt.title("Energy density")
plt.grid()
plt.tight_layout()

plt.show()
