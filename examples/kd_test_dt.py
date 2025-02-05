import numpy as np
import eigrene
from matplotlib import pyplot as plt
import math
import json
import time


config = {}

Lx = 0.1
Ly = 0.1
nx = 100
ny = 100

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
    plasma_dens[i] = 1.0e20

for i in range(ncells):
    plasma_temp[i] = 1.0

for i in range(ncells):
    plasma_vel_x[i] = 0.1

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
        "n_particles": 1000000,
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

inputstring = json.dumps(config)

start_time = time.time()
outputstring = eigrene.simulate(inputstring)
end_time = time.time()

print("Simulation took", end_time - start_time, "seconds")

output = json.loads(outputstring)

stats = output["stats"]
print(stats)

neutral_density_exact = np.array(output["density"])

neutral_density_exact_reshaped = np.reshape(neutral_density_exact, (nx, ny))

plt.figure()
c = plt.imshow(np.transpose(neutral_density_exact_reshaped), extent=[0, Lx, 0, Ly])
plt.title("Neutral particle density, pure MC")
plt.colorbar(c)
plt.xlabel("x")
plt.ylabel("y")

errors = []
dt_array = [1.0e-3, 1.0e-4, 1.0e-5, 1e-6, 1e-7]

#test with decreasing dt
for dt in dt_array:
    config["dt"] = dt

    inputstring = json.dumps(config)

    # write inputstring to file
    # with open("build/input.json", "w") as f:
    #     f.write(inputstring)

    start_time = time.time()
    outputstring = eigrene.simulate(inputstring)
    end_time = time.time()

    print("Simulation took", end_time - start_time, "seconds")

    output = json.loads(outputstring)

    stats = output["stats"]
    print(stats)

    neutral_density = np.array(output["density"])

    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

    fig, ax = plt.subplots(1,2)
    c = ax[0].imshow(np.transpose(neutral_density_reshaped), extent=[0, Lx, 0, Ly])
    ax[0].set_title("Neutral particle density, dt = " + str(dt))
    fig.colorbar(c, ax=ax[0])
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")

    diff = neutral_density_exact_reshaped - neutral_density_reshaped
    c = ax[1].imshow(np.transpose(diff/neutral_density_exact_reshaped), extent=[0, Lx, 0, Ly])
    ax[1].set_title("Difference")
    fig.colorbar(c, ax=ax[1])
    ax[1].set_xlabel("x")
    ax[1].set_ylabel("y")

    errors.append(np.linalg.norm(diff)/np.linalg.norm(neutral_density_exact_reshaped))


plt.figure()
plt.loglog(dt_array, errors, 'o-')
plt.title("Relative error vs dt")
plt.xlabel("dt")
plt.ylabel("Error")

plt.show()


