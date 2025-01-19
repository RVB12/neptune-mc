import numpy as np
import eiron2
from matplotlib import pyplot as plt
import math
import json

config = {}

Lx = 1
Ly = 1
nx = 400
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
    plasma_dens[i] = 1.0e13

for i in range(ncells):
    plasma_temp[i] = 1.0

for i in range(ncells):
    plasma_vel_x[i] = 1

config["background"] = {
    "plasma_dens": plasma_dens.tolist(),
    "plasma_temp": plasma_temp.tolist(),
    "plasma_vel_x": plasma_vel_x.tolist(),
    "plasma_vel_y": plasma_vel_y.tolist()
}


config["sources"] = [
    {
        "type": "line",
        "pos_a": [0, 0],
        "pos_b": [0, 1],
        "n_particles": 10000,
    }
]

config["max_steps"] = 1000
config["seed"] = 0

inputstring = json.dumps(config)

outputstring = eiron2.simulate(inputstring)

output = json.loads(outputstring)

neutral_density = np.array(output["density"])

x = np.linspace(0,Lx,nx)

plt.plot(x, neutral_density)
plt.title("Neutral particle density")
plt.xlabel("x")
plt.show()
