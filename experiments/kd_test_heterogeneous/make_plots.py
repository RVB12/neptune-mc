import numpy as np
from matplotlib import pyplot as plt
import math
import json
import time

Lx = 2
Ly = 1
nx = 101
ny = 1

x_pts = np.linspace(0,Lx,nx)

dt_array = np.logspace(-7, -2, 8, endpoint=True)

f = open("outputs/mc_ref.json", "r")
output = json.load(f)
neutral_density_exact = np.array(output["density"])

errors = []
runtimes = []


for dt in dt_array:

    f = open("outputs/kd_dt_" + str(dt) + ".json", "r")
    output = json.load(f)
    neutral_density = np.array(output["density"])

    plt.figure()
    plt.plot(x_pts, neutral_density_exact, label="MC")
    plt.plot(x_pts, neutral_density, label="KD")
    plt.title("dt = " + str(dt))
    plt.xlabel("x")
    plt.ylabel("density")
    plt.legend()

    plt.savefig("outputs/kd_dt_" + str(dt) + ".png")
    plt.close()

    error = np.linalg.norm(neutral_density - neutral_density_exact, 1)/np.linalg.norm(neutral_density_exact, 1)

    errors.append(error)

    runtimes.append(output["stats"]["total_time"])


plt.figure()
plt.loglog(dt_array, errors)
plt.title("Error vs dt")
plt.xlabel("dt")
plt.ylabel("Error")
plt.savefig("outputs/error_vs_dt.png")
plt.close()

# Plot runtime
plt.figure()
plt.loglog(dt_array, runtimes, 'o-')
plt.title("Runtime vs dt")
plt.xlabel("dt")
plt.ylabel("Runtime")

plt.savefig("outputs/runtime.png")
plt.close()

