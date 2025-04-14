import numpy as np
from matplotlib import pyplot as plt
import math
import json
import time

Lx = 20
Ly = 20
Lz = 20
nx = 201
ny = 1
nz = 1

x_pts = np.linspace(0,Lx,nx)

f = open("outputs/config0.1.json", "r")
output = json.load(f)

neutral_density_exact = np.array(output["density"])

eps_array = [0.1, 1, 10]

for eps in eps_array:
    f = open("outputs/config" + str(eps) + "_ref.json", "r")
    output = json.load(f)
    neutral_density_exact = np.array(output["density"])

    f = open("outputs/config" + str(eps) + ".json", "r")
    output = json.load(f)
    neutral_density = np.array(output["density"])

    plt.figure()
    plt.plot(x_pts, neutral_density_exact, label="MC")
    plt.plot(x_pts, neutral_density, label="KD")
    plt.title("epsilon = " + str(eps))
    plt.xlabel("x")
    plt.ylabel("density")
    plt.legend()

    plt.savefig("outputs/mc_eps_" + str(eps) + ".png")
    plt.close()