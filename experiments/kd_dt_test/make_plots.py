import numpy as np
from matplotlib import pyplot as plt
import math
import json
import time

nx = 21
ny = 21
nz = 1
Lx = 1
Ly = 1
Lz = 1

f = open("outputs/reference.json", "r")
output = json.load(f)

neutral_density_exact = np.array(output["density"])

neutral_density_exact_reshaped = np.reshape(neutral_density_exact, (nx, ny))

runtime_exact = output["stats"]["total_time"]

plt.figure()
c = plt.imshow(np.transpose(neutral_density_exact_reshaped), extent=[0, Lx, 0, Ly])
plt.title("Neutral particle density, pure MC")
plt.colorbar(c)
plt.xlabel("x")
plt.ylabel("y")

plt.savefig("outputs/neutral_density_exact.png")
plt.close()

errors = []
runtime = []
dt_array = np.logspace(-3, -7, num = 11, endpoint=True, base=10.0)

for dt in dt_array:

    f = open("outputs/kd_dt_" + str(dt) + ".json", "r")
    output = json.load(f)

    runtime.append(output["stats"]["total_time"])

    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))


    plt.figure()
    c = plt.imshow(np.transpose(neutral_density_reshaped), extent=[0, Lx, 0, Ly])
    plt.title("Neutral particle density, KD, dt = " + str(dt))
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("outputs/neutral_density_KD_dt_" + str(dt) + ".png")
    plt.close()

    diff = neutral_density_reshaped - neutral_density_exact_reshaped
    errors.append(np.linalg.norm(diff)/np.linalg.norm(neutral_density_exact_reshaped))

    plt.figure()
    c = plt.imshow(np.abs(np.transpose(diff))/np.max(neutral_density_exact_reshaped), extent=[0, Lx, 0, Ly])
    plt.title("Neutral particle density relative error")
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("outputs/err_KD_dt_" + str(dt) + ".png")
    plt.close()




# Plot errors
plt.figure()
plt.loglog(dt_array, errors, 'o-')
plt.title("Relative error vs dt")
plt.xlabel("dt")
plt.ylabel("Error")

plt.savefig("outputs/errors.png")
plt.close()

# Plot runtime
plt.figure()
plt.loglog(dt_array, runtime, 'o-')
plt.loglog(dt_array, runtime_exact + 0*dt_array)
plt.title("Runtime vs dt")
plt.xlabel("dt")
plt.ylabel("Runtime")

plt.savefig("outputs/runtime.png")
plt.close()

