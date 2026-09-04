# ---
# jupyter:
#   jupytext:
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

# %%
# # %matplotlib qt
# %matplotlib inline
import numpy as np
from matplotlib import pyplot as plt
import h5py
import math
import subprocess
import time
import os
import matplotlib.tri as tri

from fort_files_io import read_vertices, read_triangles, read_background, convert_background_to_neutral_grid

base_dir = "Slab_Case_Only_Atoms/Kinetic/run_ao"

vertices = read_vertices(base_dir + "/fort.33")
triangles = read_triangles(base_dir + "/fort.34")
background = read_background(base_dir + "/b2fstati")
background = convert_background_to_neutral_grid(background, base_dir + "/fort.35")

# Create the triangulation
x = np.array([v[0] for v in vertices])
y = np.array([v[1] for v in vertices])

ncells = len(triangles)
xy = np.concatenate((x[:, None], y[:, None]), axis=1)
triangulation = tri.Triangulation(xy[:, 0], xy[:, 1], triangles)

n_particles = 10000


# %%
# Plot the triangulation
plt.figure(figsize=(8, 8))
plt.triplot(triangulation, color='gray', lw=0.5)
plt.title('Triangulation')

# Plot the background for each species
for species, data in background.items():
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.tripcolor(triangulation, data["density"], norm="log")
    plt.colorbar(label=f'{species} Density')
    plt.title(f'{species} Density')
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.subplot(1, 2, 2)
    plt.tripcolor(triangulation, data["temp"])
    plt.colorbar(label=f'Temperature (eV)')
    plt.title(f'T {species}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.tight_layout()

# %%
# Make Neptune input file
import json
amjuel_data = json.load(open("../databases/amjuel.json", "r"))

with h5py.File(os.path.join("input.h5"), "w") as f:

    # Background
    bg = f.create_group("background")
    Dplus = bg.create_group("D+")
    Dplus["density"] = background["D+"]["density"]
    Dplus["density"].attrs["units"] = "m^-3"
    Dplus["velocity_x"] = background["D+"]["vel_x"]
    Dplus["velocity_x"].attrs["units"] = "m/s"
    Dplus["velocity_y"] = background["D+"]["vel_y"]
    Dplus["velocity_y"].attrs["units"] = "m/s"
    Dplus["velocity_z"] = background["D+"]["vel_z"]
    Dplus["velocity_z"].attrs["units"] = "m/s"
    Dplus["temperature"] = background["D+"]["temp"]
    Dplus["temperature"].attrs["units"] = "eV"
    electron = bg.create_group("e")
    electron["density"] = background["e"]["density"]
    electron["density"].attrs["units"] = "m^-3"
    electron["temperature"] = background["e"]["temp"]
    electron["temperature"].attrs["units"] = "eV"
    electron["velocity_x"] = np.zeros_like(background["D+"]["vel_x"])
    electron["velocity_y"] = np.zeros_like(background["D+"]["vel_y"])
    electron["velocity_z"] = np.zeros_like(background["D+"]["vel_z"])
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
    coll0["model"] = "amjuel"
    coll0["coef_matrix"] = np.array(amjuel_data["lookup"]["H.3"]["3.1.8"][0]["data_blocks"][0]["matrix_t_by_e"]).flatten()
    coll1 = collisions.create_group("1")
    coll1["background_species"] = "D+"
    coll1["particle_species"] = "D"
    coll1["type"] = "ionisation"
    coll1["model"] = "amjuel"
    coll1["coef_matrix"] = np.array(amjuel_data["lookup"]["H.4"]["2.1.5"][0]["data_blocks"][0]["matrix_t_by_e"]).flatten()
    # Sources
    sources = f.create_group("sources")
    sources["n_particles"] = n_particles
    sources["stratification"] = "source_strength"
    # Some volume recombination sources
    source = sources.create_group("R_vol")
    source["type"] = "amjuel_recombination"
    source["background_species"] = "D+"
    source["particle_species"] = "D"
    source["cell_indices"] = np.arange(ncells, dtype=np.int64)
    source["coef_matrix"] = np.array(amjuel_data["lookup"]["H.4"]["2.1.8"][0]["data_blocks"][0]["matrix_t_by_e"]).flatten()
    # Estimators
    estimators = f.create_group("estimators")
    density_estimator = estimators.create_group("density")
    density_estimator["species"] = "D"
    density_estimator["qoi"] = "density"
    density_estimator["estimator_type"] = "track_length"
    density_estimator["estimate_variance"] = 0
    velocity_x_estimator = estimators.create_group("velocity_x")
    velocity_x_estimator["species"] = "D"
    velocity_x_estimator["qoi"] = "velocity_x"
    velocity_x_estimator["estimator_type"] = "track_length"
    velocity_x_estimator["estimate_variance"] = 0
    velocity_y_estimator = estimators.create_group("velocity_y")
    velocity_y_estimator["species"] = "D"
    velocity_y_estimator["qoi"] = "velocity_y"
    velocity_y_estimator["estimator_type"] = "track_length"
    velocity_y_estimator["estimate_variance"] = 0
    velocity_z_estimator = estimators.create_group("velocity_z")
    velocity_z_estimator["species"] = "D"
    velocity_z_estimator["qoi"] = "velocity_z"
    velocity_z_estimator["estimator_type"] = "track_length"
    velocity_z_estimator["estimate_variance"] = 0
    energy_estimator = estimators.create_group("energy")
    energy_estimator["species"] = "D"
    energy_estimator["qoi"] = "energy"
    energy_estimator["estimator_type"] = "track_length"
    energy_estimator["estimate_variance"] = 0
    # Boundary conditions (specular everywhere)
    boundary_conditions = f.create_group("boundary_conditions")
    bc = boundary_conditions.create_group(f"bc")
    # bc["type"] = "specular"
    bc["type"] = "trim" 
    bc["trim_filename"] = "/home/thijs/Documents/Git/Monte-Carlo/neptune/databases/TRIM/D_on_W"
    # Simulation parameters
    params = f.create_group("params")
    params["max_steps"] = 1000000
    params["n_batches"] = 1

# %%
# Run the simulation
start_time = time.time()
output = subprocess.check_output(["neptune-mc"]).decode("utf-8")
end_time = time.time()
print(output)
print(f"Simulation completed in {end_time - start_time:.2f} seconds")

# %%
e = 1.602176634 * 1e-19 # Elementary charge
Da_to_kg = 1.66053906892 * 1e-27 # Dalton (unified atomic mass unit) to kg
m = 2. * Da_to_kg # Mass of deuterium
eV_to_J = 1.602176634 * 1e-19 # Electronvolt to Joules

# Now make some plots of the results
with h5py.File("output.h5", "r") as f:

    D_density = f["D"]["track_length"]["density"][:]
    D_velocity_x = f["D"]["track_length"]["velocity_x"][:]
    D_velocity_y = f["D"]["track_length"]["velocity_y"][:]
    D_velocity_z = f["D"]["track_length"]["velocity_z"][:]

    D_energy = f["D"]["track_length"]["energy"][:]

    for index, value in np.ndenumerate(D_density):
        if value < 1e-10:
            D_density[index] = 1e-10
            D_velocity_x[index] = 0.0
            D_velocity_y[index] = 0.0
            D_velocity_z[index] = 0.0
            D_energy[index] = 0.0


    D_temperature = 2.0/3.0 * m / e * (D_energy/D_density - 0.5 * ((D_velocity_x/D_density )**2 + (D_velocity_y/D_density )**2+ (D_velocity_z/D_density )**2 ))

    plt.figure()
    plt.tripcolor(triangulation,D_density, shading='flat', cmap="viridis", norm="log")
    plt.colorbar(label="Density (m^-3)")
    plt.title("D Density")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,D_temperature, shading='flat', cmap="viridis")
    plt.colorbar(label="Temperature (eV)")
    plt.title("D Temperature")
    plt.xlabel("r")
    plt.ylabel("z")

