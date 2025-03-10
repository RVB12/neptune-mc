import numpy as np
from matplotlib import pyplot as plt
import math
import json
import time

Lx = 20
Ly = 20
Lz = 20
nx = 101
ny = 1
nz = 1

x_pts = np.linspace(0,Lx,nx)

f = open("outputs/reference.json", "r")
output = json.load(f)

neutral_density_exact = np.array(output["density"])

runtime_exact = output["stats"]["total_time"]

plt.figure()
plt.plot(x_pts, neutral_density_exact)
plt.title("Neutral particle density, pure MC")
plt.xlabel("x")
plt.ylabel("density")

plt.savefig("outputs/neutral_density_exact.png")
plt.close()

errors = []
runtime = []
dt_array = np.array([1.0, 0.1, 0.01, 0.001])

for dt in dt_array:

    f = open("outputs/kd_dt_" + str(dt) + ".json", "r")
    output = json.load(f)

    runtime.append(output["stats"]["total_time"])

    neutral_density = np.array(output["density"])

    plt.figure()
    plt.plot(x_pts, neutral_density, label="KD")
    plt.plot(x_pts, neutral_density_exact, label="MC")
    plt.legend()
    plt.title("Neutral particle density, KD, dt = " + str(dt))
    plt.xlabel("x")
    plt.ylabel("density")

    plt.savefig("outputs/neutral_density_KD_dt_" + str(dt) + ".png")
    plt.close()

    diff = neutral_density - neutral_density_exact
    errors.append(np.linalg.norm(diff)/np.linalg.norm(neutral_density_exact))

    plt.figure()
    plt.plot(x_pts, np.abs(diff)/np.max(neutral_density_exact))
    plt.title("Neutral particle density relative error")
    plt.xlabel("x")
    plt.ylabel("density")

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

