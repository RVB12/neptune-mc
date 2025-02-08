import numpy as np
import math
import json

config = {}

Lx = 20
Ly = 20
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


for i in range(ncells):
    plasma_dens[i] = 1.0e17

for i in range(ncells):
    plasma_temp[i] = 1.0e-8

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
        "n_particles": 100000,
        "intensity": 1.0e10,
        "pos": [10, 10],
        "vel_mean": [0, 0],
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
