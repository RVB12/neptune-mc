import numpy as np
import math
import json


config = {}

Lx = 2
Ly = 1
Lz = 1
nx = 101
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
        "n_particles": 10000,
        "intensity": 1.0e10,
        "pos": [0.9, 0.5, 0.5],
        "vel_mean": [-1, 0, 0],
        "vel_stddev": 1000,
    }
]

config["boundary_conditions"] = [
    "absorbing",
    "reflective",
    "reflective",
    "reflective",
    "reflective",
    "reflective"
]

config["t_final"] = 0.1

config["max_steps"] = 10000000
config["seed"] = 0

m_p = 1.67e-27
eV_to_Joules = 1.602176634e-19

T = 1.0e8 * m_p / eV_to_Joules
dens_left = 1.0e8 / (3.2e-15 * np.sqrt(T/0.026))
dens_right = 1.0e6 / (3.2e-15 * np.sqrt(T/0.026))

plasma_dens = np.logspace(np.log10(dens_left), np.log10(dens_right), nx)

for i in range(ncells):
    plasma_temp[i] = T

for i in range(ncells):
    plasma_vel_x[i] = 500

for i in range(ncells):
    plasma_vel_y[i] = 0

config["background"] = {
    "plasma_dens": plasma_dens.tolist(),
    "plasma_temp": plasma_temp.tolist(),
    "plasma_vel_x": plasma_vel_x.tolist(),
    "plasma_vel_y": plasma_vel_y.tolist(),
    "plasma_vel_z": plasma_vel_z.tolist()
}


inputstring = json.dumps(config)
with open("inputs/mc_ref.json", "w") as f:
    f.write(inputstring)

dt_array = np.logspace(-7, -2, 8, endpoint=True)

for dt in dt_array:
    config["dt"] = dt
    inputstring = json.dumps(config)
    with open("inputs/kd_dt_" + str(dt) + ".json", "w") as f:
        f.write(inputstring)


