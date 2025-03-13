import numpy as np
import eigrene
from matplotlib import pyplot as plt
import math
import json

config = {}

Lx = 0.1
Ly = 0.1
Lz = 0.1
nx = 100
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

for i in range(ncells):
    plasma_dens[i] = 1.0e20

for i in range(ncells):
    plasma_temp[i] = 1.0

for i in range(ncells):
    plasma_vel_x[i] = 10000

config["background"] = {
    "plasma_dens": plasma_dens.tolist(),
    "plasma_temp": plasma_temp.tolist(),
    "plasma_vel_x": plasma_vel_x.tolist(),
    "plasma_vel_y": plasma_vel_y.tolist(),
    "plasma_vel_z": plasma_vel_z.tolist()
}


config["sources"] = [
    {
        "type": "point",
        "intensity": 5.0e21,
        "pos": [0.05, 0.05, 0.05],
        "vel_mean": [0, 0, 0],
        "vel_stddev": 10000,
        "n_particles": 100000,
    }
]

config["t_final"] = 1.0e-4

config["boundary_conditions"] = [
    "reflective",
    "reflective",
    "reflective",
    "reflective",
    "reflective",
    "reflective"
]

config["max_steps"] = 1000
config["seed"] = 0

inputstring = json.dumps(config)

outputstring = eigrene.simulate(inputstring)

output = json.loads(outputstring)

neutral_density = np.array(output["density"])

x = np.linspace(0,Lx,nx)

plt.plot(x, neutral_density)
plt.title("Neutral particle density")
plt.xlabel("x")
plt.show()
