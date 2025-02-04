import numpy as np
import eiron2
from matplotlib import pyplot as plt
import math
import json
import time


config = {}

Lx = 0.1
Ly = 0.1
nx = 10
ny = 10

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
    plasma_dens[i] = 1.0e21

for i in range(ncells):
    plasma_temp[i] = 1.0e8

for i in range(ncells):
    plasma_vel_x[i] = 1.0e1

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
        "type": "line",
        "pos_a": [0, 0.025],
        "pos_b": [0, 0.075],
        "n_particles": 10000,
        "intensity": 1.0e10,
        "angle_min": -math.pi/2,
        "angle_max": math.pi/2,
    }
]

config["boundary_conditions"] = [
    "reflective",
    "reflective",
    "absorbing",
    "reflective"
]

config["max_steps"] = 10000000
config["seed"] = 0

#test with increasing collision rate
for background_density in [1.e16, 1.0e17, 1.0e18, 1.0e19, 1.0e20]:
    config["dt"] = 1.0e-3

    for i in range(ncells):
        plasma_dens[i] = background_density

    config["background"]["plasma_dens"] = plasma_dens.tolist()

    inputstring = json.dumps(config)

    # # write inputstring to file
    # with open("build/input.json", "w") as f:
    #     f.write(inputstring)

    outputstring = eiron2.simulate(inputstring)

    output = json.loads(outputstring)

    stats = output["stats"]
    print(stats)
    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

    plt.figure()
    c = plt.imshow(np.transpose(neutral_density_reshaped), extent=[0, Lx, 0, Ly])
    plt.title("Neutral particle density kd, background = " + str(background_density))
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")

    del config["dt"]
    inputstring = json.dumps(config)
    outputstring = eiron2.simulate(inputstring)
    output = json.loads(outputstring)

    stats = output["stats"]
    print(stats)
    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

    plt.figure()
    c = plt.imshow(np.transpose(neutral_density_reshaped), extent=[0, Lx, 0, Ly])
    plt.title("Neutral particle density, background = " + str(background_density))
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")

plt.show()


