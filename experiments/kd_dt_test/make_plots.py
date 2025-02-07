import numpy as np
import eigrene
from matplotlib import pyplot as plt
import math
import json
import time

nx = 21
ny = 21
Lx = 1
Ly = 1

f = open("outputs/reference.json", "r")
output = json.load(f)

neutral_density_exact = np.array(output["density"])

neutral_density_exact_reshaped = np.reshape(neutral_density_exact, (nx, ny))

plt.figure()
c = plt.imshow(np.transpose(neutral_density_exact_reshaped), extent=[0, Lx, 0, Ly])
plt.title("Neutral particle density, pure MC")
plt.colorbar(c)
plt.xlabel("x")
plt.ylabel("y")

plt.savefig("outputs/neutral_density_exact.png")

errors = []
runtime = []
dt_array = [1.0e-3, 1.0e-4, 1.0e-5, 1e-6, 1.0e-7, 1.0e-8]

for dt in dt_array:

    f = open("outputs/kd_dt_" + str(dt) + ".json", "r")
    output = json.load(f)

    runtime.append(output["stats"]["total_time"])

    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))


    plt.figure()
    c = plt.imshow(np.transpose(neutral_density_reshaped), extent=[0, Lx, 0, Ly])
    plt.title("Neutral particle density, pure MC")
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("outputs/neutral_density_KD_dt_" + str(dt) + ".png")

    diff = neutral_density_reshaped - neutral_density_exact_reshaped
    errors.append(np.linalg.norm(diff)/np.linalg.norm(neutral_density_exact_reshaped))

    plt.figure()
    c = plt.imshow(np.abs(np.transpose(diff))/np.max(neutral_density_exact_reshaped), extent=[0, Lx, 0, Ly])
    plt.title("Neutral particle density, pure MC")
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("outputs/err_KD_dt_" + str(dt) + ".png")




# Plot errors
plt.figure()
plt.loglog(dt_array, errors, 'o-')
plt.title("Relative error vs dt")
plt.xlabel("dt")
plt.ylabel("Error")

plt.savefig("outputs/errors.png")

# Plot runtime
plt.figure()
plt.loglog(dt_array, runtime, 'o-')
plt.title("Runtime vs dt")
plt.xlabel("dt")
plt.ylabel("Runtime")

plt.savefig("outputs/runtime.png")

