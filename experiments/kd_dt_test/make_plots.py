import numpy as np
import eigrene
from matplotlib import pyplot as plt
import math
import json
import time

nx = 100
ny = 100
Lx = 0.1
Ly = 0.1

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
dt_array = [1.0e-3, 1.0e-4, 1.0e-5, 1e-6, 1.0e-7]

for dt in dt_array:

    f = open("outputs/kd_dt_" + str(dt) + ".json", "r")
    output = json.load(f)

    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

    plt.figure()
    c = plt.imshow(np.transpose(neutral_density_reshaped), extent=[0, Lx, 0, Ly])
    plt.title("Neutral particle density, KD, dt = " + str(dt))
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("outputs/neutral_density_KD_dt_" + str(dt) + ".png")

    diff = neutral_density - neutral_density_exact
    errors.append(np.linalg.norm(diff)/np.linalg.norm(neutral_density_exact_reshaped))


plt.figure()
plt.loglog(dt_array, errors, 'o-')
plt.title("Relative error vs dt")
plt.xlabel("dt")
plt.ylabel("Error")

plt.savefig("outputs/errors.png")