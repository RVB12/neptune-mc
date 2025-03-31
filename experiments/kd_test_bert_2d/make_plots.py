import numpy as np
from matplotlib import pyplot as plt
import math
import json
import time

Lx = 40
Ly = 40
nx = 201
ny = 201

eps_array = [0.1, 1, 10]

for eps in eps_array:
    f = open("outputs/config" + str(eps) + "_ref.json", "r")
    output = json.load(f)
    neutral_density_exact = np.array(output["density"])
    neutral_density_exact_reshaped = np.reshape(neutral_density_exact, (nx, ny))

    f = open("outputs/config" + str(eps) + ".json", "r")
    output = json.load(f)
    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    # x = np.linspace(0, Lx, nx)
    # y = np.linspace(0, Ly, ny)

    # X, Y = np.meshgrid(x, y)

    # ax.plot_surface(X, Y, neutral_density_exact_reshaped, cmap='viridis')
    # ax.set_title("Neutral particle density, pure MC, eps = " + str(eps))

    fig, ax = plt.subplots(2,2)
    c = ax[0,0].imshow(neutral_density_exact_reshaped, extent=[0, Lx, 0, Ly])
    ax[0,0].set_title("Neutral particle density, pure MC")
    ax[0,0].set_xlabel("x")
    ax[0,0].set_ylabel("y")
    fig.colorbar(c, ax=ax[0,0])

    c = ax[0,1].imshow(neutral_density_reshaped, extent=[0, Lx, 0, Ly])
    ax[0,1].set_title("Neutral particle density, KD")
    ax[0,1].set_xlabel("x")
    ax[0,1].set_ylabel("y")
    fig.colorbar(c, ax=ax[0,1])

    diff = neutral_density_reshaped - neutral_density_exact_reshaped
    c = ax[1,0].imshow(diff, extent=[0, Lx, 0, Ly])
    ax[1,0].set_title("Difference")
    ax[1,0].set_xlabel("x")
    ax[1,0].set_ylabel("y")
    fig.colorbar(c, ax=ax[1,0])


    err = np.linalg.norm(diff) / np.linalg.norm(neutral_density_exact)
    print("Error for epsilon = ", eps, " is ", err)

    plt.savefig("outputs/mc_eps_" + str(eps) + ".png")

plt.show()