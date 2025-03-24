import numpy as np
import math
import json


config = {}

Lx = 35
Ly = 40
Lz = 40
nx = 51
ny = 1
nz = 1

config["mesh"] = {
    "type": "voxel",
    "Lx": Lx,
    "Ly": Ly,
    "Lz": Lz,
    "nx": nx,
    "ny": ny,
    "nz": nz
}

ncells = nx * ny * nz

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)
plasma_vel_z = np.zeros(ncells, dtype=np.float64)

config["sources"] = [
    {
        "type": "point",
        "n_particles": 5000000,
        "intensity": 1.0e10,
        "pos": [15, 20, 20],
        "vel_mean": [-5, 0, 0],
        "vel_stddev": 1,
    },
    {
        "type": "point",
        "n_particles": 5000000,
        "intensity": 1.0e10,
        "pos": [15, 20, 20],
        "vel_mean": [15, 0, 0],
        "vel_stddev": 1,
    }
]

config["boundary_conditions"] = [
    "reflective",
    "reflective",
    "reflective",
    "reflective",
    "reflective",
    "reflective"
]

config["t_final"] = 1.0
config["dt"] = 1.0

config["max_steps"] = 10000000
config["seed"] = 0


eps_array = [0.1, 1, 10]
m_p = 1.67e-27
eV_to_Joules = 1.602176634e-19

for epsilon in eps_array:
    T = 1/(epsilon*epsilon) * m_p / eV_to_Joules
    dens = 1/(epsilon*epsilon) / (3.2e-15 * np.sqrt(T/0.026))

    for i in range(ncells):
        plasma_dens[i] = dens

    for i in range(ncells):
        plasma_temp[i] = T

    for i in range(ncells):
        plasma_vel_x[i] = 0

    for i in range(ncells):
        plasma_vel_y[i] = 0

    config["background"] = {
        "plasma_dens": plasma_dens.tolist(),
        "plasma_temp": plasma_temp.tolist(),
        "plasma_vel_x": plasma_vel_x.tolist(),
        "plasma_vel_y": plasma_vel_y.tolist(),
        "plasma_vel_z": plasma_vel_z.tolist()
    }

    config["dt"] = 1.0
    inputstring = json.dumps(config)
    with open("inputs/config" + str(epsilon) + ".json", "w") as f:
        f.write(inputstring)
    del config["dt"]

    inputstring = json.dumps(config)
    with open("inputs/config" + str(epsilon) + "_ref.json", "w") as f:
        f.write(inputstring)


