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
#     display_name: .venv (3.12.3)
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

from fort_files_io import read_eirene_output, read_neighbor_list, read_vertices, read_triangles, read_background, convert_background_to_neutral_grid, read_magnetic_field

# base_dir = "Slab_Case_Only_Atoms/Kinetic/baserun"
# run_dir = "Slab_Case_Only_Atoms/Kinetic/run_ao"

base_dir = "CMOD-data"
run_dir = "CMOD-data"

vertices = read_vertices(run_dir + "/fort.33")
triangles = read_triangles(run_dir + "/fort.34")
neighbors, neighbors_bc = read_neighbor_list(run_dir + "/fort.35")
background = read_background(run_dir + "/b2fstati")
magnetic_field = read_magnetic_field(base_dir + "/b2fgmtry")
background = convert_background_to_neutral_grid(background, magnetic_field, run_dir + "/fort.35")
eirene_output = read_eirene_output(run_dir + "/fort.46", run_dir + "/fort.35")

# default neighbors_bc to 6 if the neighbor is -1 (i.e. no neighbor)
# this overrides the neptune behavior of defaulting to 0, which would be an absorbing bc
for i in range(len(triangles)):
    for j in range(3):
        if neighbors[i][j] == -1 and neighbors_bc[i][j] == -1:
            neighbors_bc[i][j] = 6

# Create the triangulation
x = np.array([v[0] for v in vertices])
y = np.array([v[1] for v in vertices])

ncells = len(triangles)
xy = np.concatenate((x[:, None], y[:, None]), axis=1)
triangulation = tri.Triangulation(xy[:, 0], xy[:, 1], triangles)

n_particles = 100000


# %%
# %matplotlib inline
# Plot the triangulation
plt.figure(figsize=(8, 8))
plt.triplot(triangulation, color='gray', lw=0.5)
plt.title('Triangulation')
plt.xlabel('r')
plt.ylabel('z')


# Plot the background for each species
for species, data in background.items():
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.tripcolor(triangulation, data["density"], norm="log")
    plt.colorbar(label=f'{species} Density')
    plt.title(f'{species} Density')
    plt.xlabel('r')
    plt.ylabel('z')

    plt.subplot(1, 2, 2)
    plt.tripcolor(triangulation, data["temp"], norm="log")
    plt.colorbar(label=f'Temperature (eV)')
    plt.title(f'T {species}')
    plt.xlabel('r')
    plt.ylabel('z')
    plt.tight_layout()

    if species == "D+":
        plt.figure()
        plt.tripcolor(triangulation, data["vel_x"])
        plt.colorbar(label=f'Velocity X (m/s)')
        plt.title(f'Velocity X (radial) {species}')
        plt.xlabel('r')
        plt.ylabel('z')

        plt.figure()
        plt.tripcolor(triangulation, data["vel_y"])
        plt.colorbar(label=f'Velocity Y (m/s)')
        plt.title(f'Velocity Y (toroidal) {species}')
        plt.xlabel('r')
        plt.ylabel('z')

        plt.figure()
        plt.tripcolor(triangulation, data["vel_z"])
        plt.colorbar(label=f'Velocity Z (m/s)')
        plt.title(f'Velocity Z (poloidal) {species}')
        plt.xlabel('r')
        plt.ylabel('z')

        # Quiver plot for velocity vectors

        cell_centers_x = np.mean(x[triangles], axis=1)
        cell_centers_y = np.mean(y[triangles], axis=1)
        plt.figure()
        plt.quiver(cell_centers_x, cell_centers_y, data["vel_x"], data["vel_z"], scale=1e5)
        


        

# %%
# # %matplotlib qt
# # Show the boundary conditions with different colors

# bc_colors = np.max(neighbors_bc, axis=1)
# bc_plot = np.ma.masked_where(bc_colors == 0, bc_colors)
# plt.figure(figsize=(8, 8))
# plt.triplot(triangulation, color='gray', lw=0.8)
# cmap = plt.cm.get_cmap('tab10').copy()
# cmap.set_bad((1, 1, 1, 0))  # transparent where bc_colors == -1
# plt.tripcolor(triangulation, facecolors=bc_plot, cmap=cmap, alpha=0.7, lw=0.1, vmin=0, vmax=10)
# plt.title('Bound0ary Conditions')
# plt.xlabel('r')
# plt.ylabel('z')
# plt.colorbar(label='Boundary Condition Type')
# # Also put a tiny number in the middle of each triangle to indicate the cell number
# for i, tri_indices in enumerate(triangles):
#     # only label the triangles that are not masked (i.e., have a boundary condition)
#     if bc_colors[i] == -1:
#         continue
#     tri_center_x = np.mean(x[tri_indices])
#     tri_center_y = np.mean(y[tri_indices])
#     plt.text(tri_center_x, tri_center_y, str(bc_colors[i]), fontsize=6, ha='center', va='center')

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
    mesh["type"] = "rotated_triangle"
    mesh["vertices"] = vertices
    mesh["triangles"] = np.array(triangles, dtype=np.int64)
    mesh["neighbors"] = neighbors.flatten()
    mesh["boundary_conditions"] = neighbors_bc.flatten()
    # Collisions
    collisions = f.create_group("collisions")
    coll0 = collisions.create_group("0")
    coll0["background_species"] = "D+"
    coll0["particle_species"] = "D"
    coll0["type"] = "charge_exchange"
    # coll0["model"] = "simple"
    coll0["model"] = "amjuel"
    coll0["coef_matrix"] = np.array(amjuel_data["lookup"]["H.3"]["3.1.8"][0]["data_blocks"][0]["matrix_t_by_e"]).flatten()
    coll0["output_collision_info"] = 1
    coll1 = collisions.create_group("1")
    coll1["background_species"] = "D+"
    coll1["particle_species"] = "D"
    coll1["type"] = "ionisation"
    coll1["model"] = "simple"
    # coll1["model"] = "amjuel"
    # coll1["coef_matrix"] = np.array(amjuel_data["lookup"]["H.4"]["2.1.5"][0]["data_blocks"][0]["matrix_t_by_e"]).flatten()
    coll1["output_collision_info"] = 1
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
    # Wall recycling source
    # This is a bit of a manual hack, we apply surface recycling only on the cells with boundary condition -3
    cell_indices_wall = []
    surface_indices_wall = []
    for i, bc in enumerate(neighbors_bc):
        if bc[0] == 1:
            cell_indices_wall.append(i)
            surface_indices_wall.append(0)
        if bc[1] == 1:
            cell_indices_wall.append(i)
            surface_indices_wall.append(1)
        if bc[2] == 1:
            cell_indices_wall.append(i)
            surface_indices_wall.append(2)
    cell_indices_wall = np.array(cell_indices_wall, dtype=np.int64)
    surface_indices_wall = np.array(surface_indices_wall, dtype=np.int64)
    wall_source = sources.create_group("R_wall")
    wall_source["type"] = "wall_recycling"
    wall_source["background_species"] = "D+"
    wall_source["particle_species"] = "D"
    wall_source["cell_indices"] = cell_indices_wall
    wall_source["surface_indices"] = surface_indices_wall
    wall_source["sheath_potential"] = 3.1
    wall_source["franck_condon_temp"] = 2.0
    wall_source["Rt"] = 0.9
    wall_source["Rr"] = 0.9

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
    # Boundary conditions
    # TODO: actually read the boundary conditions from the file
    boundary_conditions = f.create_group("boundary_conditions")
    # bc1: absorbing bc at the plasma core
    bc1 = boundary_conditions.create_group("bc1")
    bc1["type"] = "absorbing"
    # bc2: west wall
    bc2 = boundary_conditions.create_group("bc2")
    bc2["type"] = "trim"
    bc2["trim_filename"] = "/home/thijs/Documents/Git/Monte-Carlo/neptune/databases/TRIM/D_on_Mo"
    # bc3: east wall
    bc3 = boundary_conditions.create_group("bc3")
    bc3["type"] = "trim"
    bc3["trim_filename"] = "/home/thijs/Documents/Git/Monte-Carlo/neptune/databases/TRIM/D_on_Mo"
    # bc4, bc5, and bc6, no physical wall, just transfer to void region using transparent boundary conditions
    bc4 = boundary_conditions.create_group("bc4")
    bc4["type"] = "transparent"
    bc5 = boundary_conditions.create_group("bc5")
    bc5["type"] = "transparent"
    bc6 = boundary_conditions.create_group("bc6")
    bc6["type"] = "transparent"
    # bc7: extra boundary condition for any cell that has no neighbor
    bc7 = boundary_conditions.create_group("bc7")
    bc7["type"] = "diffuse"
    bc7["fixed_energy"] = 1.0
    bc7["absorption_coefficient"] = 0.0

    # for i in range(max(-neighbors_bc.flatten())):
    #     bc = boundary_conditions.create_group(f"bc{i}")
    #     bc["type"] = "absorbing"
        # bc["type"] = "trim" 
        # bc["trim_filename"] = "/home/thijs/Documents/Git/Monte-Carlo/neptune/databases/TRIM/D_on_W"
        # bc["type"] = "diffuse"
        # bc["fixed_energy"] = 1.0
        # bc["absorption_coefficient"] = 0.1
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
# %matplotlib inline
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

    # # Set void region to 0 density and 0 velocity
    # # This makes it easier to compare against the eirene output
    # for index, value in np.ndenumerate(background["D+"]["density"]):
    #     if value < 1e-10:
    #         D_density[index] = 1e-10
    #         D_velocity_x[index] = 0.0
    #         D_velocity_y[index] = 0.0
    #         D_velocity_z[index] = 0.0
    #         D_energy[index] = 0.0

    for index, value in np.ndenumerate(eirene_output["density"]):
        if value < 1e-10:
            eirene_output["density"][index] = 1e-10
            eirene_output["vel_x"][index] = 0.0
            eirene_output["vel_y"][index] = 0.0
            eirene_output["vel_z"][index] = 0.0
            eirene_output["energy"][index] = 0.0

    D_temperature = 2.0/3.0 * m / e * (D_energy/D_density - 0.5 * ((D_velocity_x/D_density )**2 + (D_velocity_y/D_density )**2+ (D_velocity_z/D_density )**2 ))

    plt.figure()
    plt.tripcolor(triangulation,D_density, shading='flat', cmap="viridis", norm="log", vmin=1e14, vmax=1e21)
    plt.colorbar(label="Density (m^-3)")
    plt.title("D Density (Neptune)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,eirene_output["density"], shading='flat', cmap="viridis", norm="log", vmin=1e14, vmax=1e21)
    plt.colorbar(label="Density (m^-3)")
    plt.title("D Density (Eirene)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,f["collision_info"]["charge_exchange_0"]["collision_rate"], shading='flat', cmap="viridis", norm="log")
    plt.colorbar(label="Collision Rate (s^-1)")
    plt.title("Charge Exchange Rate (Neptune)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,f["collision_info"]["ionisation_0"]["collision_rate"], shading='flat', cmap="viridis", norm="log")
    plt.colorbar(label="Collision Rate (s^-1)")
    plt.title("Ionisation Rate (Neptune)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,D_velocity_x/D_density, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity X (m/s)")
    plt.title("D Velocity X (Neptune)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,eirene_output["vel_x"]/eirene_output["density"]/(m*1000)*1e-2, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity X (m/s)")
    plt.title("D Velocity X (Eirene)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,D_velocity_y/D_density, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Y (m/s)")
    plt.title("D Velocity Y (Neptune)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,eirene_output["vel_y"]/eirene_output["density"]/(m*1000)*1e-2, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Y (m/s)")
    plt.title("D Velocity Y (Eirene)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,D_velocity_z/D_density, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Z (m/s)")
    plt.title("D Velocity Z (Neptune)")
    plt.xlabel("r")
    plt.ylabel("z")

    plt.figure()
    plt.tripcolor(triangulation,eirene_output["vel_z"]/eirene_output["density"]/(m*1000)*1e-2, shading='flat', cmap="viridis")
    plt.colorbar(label="Velocity Z (m/s)")
    plt.title("D Velocity Z (Eirene)")
    plt.xlabel("r")
    plt.ylabel("z")

