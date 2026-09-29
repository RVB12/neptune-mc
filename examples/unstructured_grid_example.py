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
#     display_name: .venv (3.12.3.final.0)
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
# Step 1, create an unstructured grid
import matplotlib.tri as tri

# First create the x and y coordinates of the points.
n_angles = 60
n_radii = 60
min_radius = 0.01
max_radius = 1
offsetr = 2.0
offsetz = 0.0

radii = np.linspace(min_radius, max_radius, n_radii)

angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
# angles[:, 1::2] += np.pi / n_angles

r = (radii * np.cos(angles)).flatten() + offsetr
z = (radii * np.sin(angles)).flatten() + offsetz

triang = tri.Triangulation(r, z)

vertices = np.concatenate((r[:, None], z[:, None]), axis=1)
triangles = triang.get_masked_triangles()
ncells = triangles.shape[0]

# Calculate the cell centers and volumes
rc = np.zeros(ncells)
zc = np.zeros(ncells)
for i in range(ncells):
    v1 = vertices[triangles[i, 0]]
    v2 = vertices[triangles[i, 1]]
    v3 = vertices[triangles[i, 2]]
    rc[i] = (v1[0] + v2[0] + v3[0]) / 3
    zc[i] = (v1[1] + v2[1] + v3[1]) / 3

# %%
# Now determine the phyical parameters on the grid

e = 1.602176634 * 1e-19 # Elementary charge
Da_to_kg = 1.66053906892 * 1e-27 # Dalton (unified atomic mass unit) to kg
m = 2. * Da_to_kg # Mass of deuterium
eV_to_J = 1.602176634 * 1e-19 # Electronvolt to Joules

# Temperature in eV, varies linearly from 10 eV at the center of the circle to 2 eV at the edge of the circle
ion_temperature = 10 - 8 * np.sqrt((rc - offsetr)**2 + (zc - offsetz)**2) / max_radius

# No drift in radial (r) or poloidal (z) direction, thermal velocity in toroidal (t) direction
vel_r = np.zeros(ncells, dtype=np.float64)
vel_z = np.zeros(ncells, dtype=np.float64)
sigma_plasma = np.sqrt(e * ion_temperature / m)
# This toroidal velocity is absurdly high, but it makes the effect much clearer
vel_t = 50000 * np.ones(ncells, dtype=np.float64)

# Uniform plasma density
ion_density = 1e20 * np.ones(ncells, dtype=np.float64)

# Plot the parameters on the grid
plt.figure()
plt.tripcolor(r,z,ion_density, triangles=triangles, shading='flat', cmap="viridis", norm="log")
plt.colorbar(label="Ion Density (m^-3)")
plt.title("Ion Density")
plt.xlabel("r")
plt.ylabel("z")

plt.figure()
plt.tripcolor(r,z,vel_t, triangles=triangles, shading='flat', cmap="viridis")
plt.colorbar(label="Ion Velocity (m/s)")
plt.title("Ion Velocity in toroidal direction")
plt.xlabel("r")
plt.ylabel("z")

plt.figure()
plt.tripcolor(r,z,ion_temperature, triangles=triangles, shading='flat', cmap="viridis")
plt.colorbar(label="Ion Temperature (eV)")
plt.title("Ion Temperature")
plt.xlabel("r")
plt.ylabel("z")

n_particles = 100000
estimate_variance = False
n_batches = 1

# %%
#
# Transform these parameters to a hdf5 file that can be read by the Monte Carlo code
#
with h5py.File(os.path.join("input.h5"), "w") as f:

    # Background
    background = f.create_group("background")
    Dplus = background.create_group("D+")
    Dplus["density"] = ion_density
    Dplus["density"].attrs["units"] = "m^-3"
    Dplus["velocity_x"] = vel_r
    Dplus["velocity_x"].attrs["units"] = "m/s"
    Dplus["velocity_y"] = vel_t
    Dplus["velocity_y"].attrs["units"] = "m/s"
    Dplus["velocity_z"] = vel_z
    Dplus["velocity_z"].attrs["units"] = "m/s"
    Dplus["temperature"] = ion_temperature
    Dplus["temperature"].attrs["units"] = "eV"
    electron = background.create_group("e")
    electron["density"] = ion_density
    electron["density"].attrs["units"] = "m^-3"
    electron["temperature"] = ion_temperature
    electron["temperature"].attrs["units"] = "eV"
    electron["velocity_x"] = np.zeros_like(vel_r)
    electron["velocity_y"] = np.zeros_like(vel_r)
    electron["velocity_z"] = np.zeros_like(vel_r)
    # Mesh
    mesh = f.create_group("mesh")
    mesh["type"] = "triangle"
    mesh["vertices"] = vertices
    mesh["triangles"] = np.array(triangles, dtype=np.int64)
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
    source["cell_indices"] = np.arange(ncells, dtype=np.int64)
    # Estimators
    estimators = f.create_group("estimators")
    density_estimator = estimators.create_group("density")
    density_estimator["species"] = "D"
    density_estimator["qoi"] = "density"
    density_estimator["estimator_type"] = "track_length"
    density_estimator["estimate_variance"] = 1 if estimate_variance else 0
    velocity_x_estimator = estimators.create_group("velocity_x")
    velocity_x_estimator["species"] = "D"
    velocity_x_estimator["qoi"] = "velocity_x"
    velocity_x_estimator["estimator_type"] = "track_length"
    velocity_x_estimator["estimate_variance"] = 1 if estimate_variance else 0
    velocity_y_estimator = estimators.create_group("velocity_y")
    velocity_y_estimator["species"] = "D"
    velocity_y_estimator["qoi"] = "velocity_y"
    velocity_y_estimator["estimator_type"] = "track_length"
    velocity_y_estimator["estimate_variance"] = 1 if estimate_variance else 0
    velocity_z_estimator = estimators.create_group("velocity_z")
    velocity_z_estimator["species"] = "D"
    velocity_z_estimator["qoi"] = "velocity_z"
    velocity_z_estimator["estimator_type"] = "track_length"
    velocity_z_estimator["estimate_variance"] = 1 if estimate_variance else 0
    energy_estimator = estimators.create_group("energy")
    energy_estimator["species"] = "D"
    energy_estimator["qoi"] = "energy"
    energy_estimator["estimator_type"] = "track_length"
    energy_estimator["estimate_variance"] = 1 if estimate_variance else 0
    energy_x_estimator = estimators.create_group("energy_x")
    energy_x_estimator["species"] = "D"
    energy_x_estimator["qoi"] = "energy_x"
    energy_x_estimator["estimator_type"] = "track_length"
    energy_x_estimator["estimate_variance"] = 1 if estimate_variance else 0
    # Boundary conditions (specular everywhere)
    boundary_conditions = f.create_group("boundary_conditions")
    bc = boundary_conditions.create_group(f"bc")
    bc["type"] = "specular"
    # Simulation parameters
    params = f.create_group("params")
    params["max_steps"] = 1000000
    params["n_batches"] = n_batches

# %%
# Run the simulation
start_time = time.time()
output = subprocess.check_output(["neptune-mc"]).decode("utf-8")
end_time = time.time()
print(output)
print(f"Simulation completed in {end_time - start_time:.2f} seconds")

# %%
# Now make some plots of the results
with h5py.File("output.h5", "r") as f:

    D_density = f["D"]["track_length"]["density"][:]
    D_velocity_x = f["D"]["track_length"]["velocity_x"][:]
    D_velocity_y = f["D"]["track_length"]["velocity_y"][:]
    D_velocity_z = f["D"]["track_length"]["velocity_z"][:]

    D_energy = f["D"]["track_length"]["energy"][:]

    D_temperature = 2.0/3.0 * m / e * (D_energy/D_density - 0.5 * ((D_velocity_x/D_density )**2 + (D_velocity_y/D_density )**2+ (D_velocity_z/D_density )**2 ))

    plt.figure()
    plt.tripcolor(r,z,D_density, triangles=triangles, shading='flat', cmap="viridis", norm="log")
    plt.colorbar(label="Density (m^-3)")
    plt.title("D Density")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(r,z,D_velocity_x/D_density, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity X (m/s)")
    plt.title("Radial Velocity")
    plt.xlabel("r")
    plt.ylabel("z")


    plt.figure()
    plt.tripcolor(r,z,D_velocity_y/D_density, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Y (m/s)")
    plt.title("Toroidal Velocity")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(r,z,D_velocity_z/D_density, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Z (m/s)")
    plt.title("Poloidal Velocity")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(r,z,D_temperature, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Temperature (eV)")
    plt.title("D Temperature")
    plt.xlabel("r")
    plt.ylabel("z")


    # Also plot arrow plot of the velocity field
    plt.figure()
    plt.tripcolor(r,z,D_density, triangles=triangles, shading='flat', cmap="viridis", norm="log")
    plt.colorbar(label="D Density (m^-3)")
    plt.title("D Density with Velocity Field")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.quiver(rc, zc, D_velocity_x/D_density, D_velocity_z/D_density, color="white", scale=3e5)

# %%
# Now rerun the simulation, but with a rotated triangle mesh
# Because of a relatively large toroidal ion velocity,
# We expect the particles to be thrown towards the right side of the domain

with h5py.File("input.h5", "a") as f:
    # Modify mesh/type to "rotated_triangle"
    f["mesh/type"][...] = "rotated_triangle"

# Rerun the simulation
start_time = time.time()
output = subprocess.check_output(["neptune-mc"]).decode("utf-8")
end_time = time.time()
print(output)
print(f"Simulation completed in {end_time - start_time:.2f} seconds")

# %%


# Now make some plots of the results
with h5py.File("output.h5", "r") as f:

    D_density = f["D"]["track_length"]["density"][:]
    D_velocity_x = f["D"]["track_length"]["velocity_x"][:]
    D_velocity_y = f["D"]["track_length"]["velocity_y"][:]
    D_velocity_z = f["D"]["track_length"]["velocity_z"][:]

    D_energy = f["D"]["track_length"]["energy"][:]

    D_temperature = 2.0/3.0 * m / e * (D_energy/D_density - 0.5 * ((D_velocity_x/D_density )**2 + (D_velocity_y/D_density )**2+ (D_velocity_z/D_density )**2 ))


    plt.figure()
    plt.tripcolor(r,z,D_density, triangles=triangles, shading='flat', cmap="viridis", norm="log")
    plt.colorbar(label="D Density (m^-3)")
    plt.title("D Density")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(r,z,D_velocity_x/D_density, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity X (m/s)")
    plt.title("Radial Velocity")
    plt.xlabel("r")
    plt.ylabel("z")


    plt.figure()
    plt.tripcolor(r,z,D_velocity_y/D_density, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Y (m/s)")
    plt.title("Toroidal Velocity")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(r,z,D_velocity_z/D_density, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Z (m/s)")
    plt.title("Poloidal Velocity")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(r,z,D_temperature, triangles=triangles, shading='flat', cmap="viridis")
    plt.colorbar(label="Temperature (eV)")
    plt.title("D Temperature")
    plt.xlabel("r")
    plt.ylabel("z")

    # Also plot arrow plot of the velocity field
    plt.figure()
    plt.tripcolor(r,z,D_density, triangles=triangles, shading='flat', cmap="viridis", norm="log")
    plt.colorbar(label="D Density (m^-3)")
    plt.title("D Density with Velocity Field")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.quiver(rc, zc, D_velocity_x/D_density, D_velocity_z/D_density, color="white", scale=6e5)

plt.show()
