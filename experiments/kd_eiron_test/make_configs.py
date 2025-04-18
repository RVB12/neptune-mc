import numpy as np
import math
import json


config = {}

Lx = 128
Ly = 128
Lz = 128
nx = 128
ny = 128
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
        "n_particles": 1000000,
        "intensity": 1.0,
        "pos": [64, 64, 64],
        "vel_stddev": 1,
        "vel_mean": [0, 0, 0],
    },
]

config["boundary_conditions"] = [
    "reflective",
    "reflective",
    "reflective",
    "reflective",
    "reflective",
    "reflective"
]

config["t_final"] = 20.0
config["dt"] = 1.0

config["max_steps"] = 10000000
config["seed"] = 0


desired_bgk_stddev = 1
desired_collision_rate = 10

def get_BGK_stddev(T):
    return np.sqrt(T * 1.602176634e-19 / 1.67e-27)

def get_collision_rate(T, dens):
    return 3.2e-15 * np.sqrt(T/0.026) * dens

for desired_collision_rate in [0.01, 0.1, 1, 10]:

    m_p = 1.67e-27
    eV_to_Joules = 1.602176634e-19

    T = desired_bgk_stddev*desired_bgk_stddev * m_p / eV_to_Joules
    dens = desired_collision_rate / (3.2e-15 * np.sqrt(T/0.026))
    
    # Calculate expected variance
    bgk_stddev = get_BGK_stddev(T)
    collision_rate = get_collision_rate(T, dens)
    print("BGK stddev: ", bgk_stddev)
    print("Collision rate: ", collision_rate)
    print("solution stddev: ", 2*bgk_stddev / np.sqrt(collision_rate))



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
    with open("inputs/config" + str(desired_collision_rate) + ".json", "w") as f:
        f.write(inputstring)
    del config["dt"]

    inputstring = json.dumps(config)
    with open("inputs/config" + str(desired_collision_rate) + "_ref.json", "w") as f:
        f.write(inputstring)


