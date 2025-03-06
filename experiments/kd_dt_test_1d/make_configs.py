import numpy as np
import math
import json

config = {}

Lx = 40
Ly = 40
nx = 101
ny = 1

config["mesh"] = {
    "type": "voxel",
    "Lx": Lx,
    "Ly": Ly,
    "nx": nx,
    "ny": ny
}

ncells = nx * ny

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)

eps_array = [0.1, 1, 10]
m_p = 1.67e-27
eV_to_Joules = 1.602176634e-19
epsilon = 1
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
    "plasma_vel_y": plasma_vel_y.tolist()
}

config["sources"] = [
    {
        "type": "point",
        "n_particles": 50000000,
        "intensity": 1.0e10,
        "pos": [20, 20],
        "vel_mean": [-15, 0],
        "vel_stddev": 1,
    },
    {
        "type": "point",
        "n_particles": 50000000,
        "intensity": 1.0e10,
        "pos": [20, 20],
        "vel_mean": [5, 0],
        "vel_stddev": 1,
    }
]

config["boundary_conditions"] = [
    "reflective",
    "reflective",
    "reflective",
    "reflective"
]

config["t_final"] = 1.0

config["max_steps"] = 10000000
config["seed"] = 0

inputstring = json.dumps(config)
with open("inputs/reference.json", "w") as f:
    f.write(inputstring)

dt_array = [1.0, 0.1, 0.01, 0.001]

#test with decreasing dt
for dt in dt_array:
    config["dt"] = dt

    inputstring = json.dumps(config)
    with open("inputs/kd_dt_" + str(dt) + ".json", "w") as f:
        f.write(inputstring)
