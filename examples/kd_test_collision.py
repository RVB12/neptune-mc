import numpy as np
import eigrene
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
    plasma_temp[i] = 1.0

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
        "n_particles": 1000,
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

errors = []
background_density_array = [1.0e16, 1.0e17, 1.0e18, 1.0e19, 1.0e20, 1.0e21, 1.0e22]

#test with increasing collision rate
for background_density in background_density_array:
    config["dt"] = 1.0e-3

    for i in range(ncells):
        plasma_dens[i] = background_density

    config["background"]["plasma_dens"] = plasma_dens.tolist()

    inputstring = json.dumps(config)

    # # write inputstring to file
    # with open("build/input.json", "w") as f:
    #     f.write(inputstring)

    outputstring = eigrene.simulate(inputstring)

    output = json.loads(outputstring)

    stats = output["stats"]
    print(stats)
    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

    fig, ax = plt.subplots(1,2)
    c = ax[0].imshow(np.transpose(neutral_density_reshaped), extent=[0, Lx, 0, Ly])
    ax[0].set_title("Neutral particle density (KD), background = " + str(background_density))
    fig.colorbar(c, ax=ax[0])
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")

    del config["dt"]
    inputstring = json.dumps(config)
    outputstring = eigrene.simulate(inputstring)
    output = json.loads(outputstring)

    stats = output["stats"]
    print(stats)
    neutral_density = np.array(output["density"])
    neutral_density_exact_reshaped = np.reshape(neutral_density, (nx, ny))

    c = ax[1].imshow(np.transpose(neutral_density_exact_reshaped), extent=[0, Lx, 0, Ly])
    ax[1].set_title("Neutral particle density, background = " + str(background_density))
    fig.colorbar(c, ax=ax[1])
    ax[1].set_xlabel("x")
    ax[1].set_ylabel("y")

    diff = neutral_density_exact_reshaped - neutral_density_reshaped
    errors.append(np.linalg.norm(diff) / np.linalg.norm(neutral_density_exact_reshaped))
    print("Error: ", errors[-1])


plt.figure()
plt.loglog(background_density_array, errors)
plt.title("Error vs background density")
plt.xlabel("Background density")
plt.ylabel("Relative error")

plt.show()


