import numpy as np
from matplotlib import pyplot as plt
import math
import json
import time

Lx = 128
Ly = 128
Lz = 128
nx = 128
ny = 128
nz = 1

for desired_collision_rate in [0.01, 0.1, 1, 10]:

    f = open("outputs/config" + str(desired_collision_rate) + "_ref.json", "r")
    output = json.load(f)
    neutral_density_exact = np.array(output["density"])
    neutral_density_exact_reshaped = np.reshape(neutral_density_exact, (nx, ny))

    f = open("outputs/config" + str(desired_collision_rate) + ".json", "r")
    output = json.load(f)
    neutral_density = np.array(output["density"])
    neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

    # make a 3D plot of the neutral density
    X, Y = np.meshgrid(np.linspace(0, Lx, nx), np.linspace(0, Ly, ny))
    fig, ax = plt.subplots(2,2, subplot_kw={"projection": "3d"})
    ax[0,0].plot_surface(X, Y, neutral_density_exact_reshaped, cmap='viridis')
    ax[0,0].set_title("Neutral particle density, pure MC")
    ax[0,0].set_xlabel("x")
    ax[0,0].set_ylabel("y")
    ax[0,0].set_zlabel("density")

    ax[0,1].plot_surface(X, Y, neutral_density_reshaped, cmap='viridis')
    ax[0,1].set_title("Neutral particle density, KD")
    ax[0,1].set_xlabel("x")
    ax[0,1].set_ylabel("y")
    ax[0,1].set_zlabel("density")
    diff = neutral_density_reshaped - neutral_density_exact_reshaped
    ax[1,0].plot_surface(X, Y, diff, cmap='viridis')
    ax[1,0].set_title("Difference")
    ax[1,0].set_xlabel("x")
    ax[1,0].set_ylabel("y")
    ax[1,0].set_zlabel("density")
    
    # expected density: normal around 64,64 with stddev = sqrt(2/desired_collision_rate * 20)

    expected_stddev = np.sqrt(1/desired_collision_rate * 20)
    expected_density = (1/Lz) * 1/(2*math.pi*expected_stddev**2) * np.exp(-((X-64)**2 + (Y-64)**2)/(2*expected_stddev**2))
    ax[1,1].plot_surface(X, Y, expected_density, cmap='viridis')
    ax[1,1].set_title("Expected density")
    ax[1,1].set_xlabel("x")
    ax[1,1].set_ylabel("y")
    ax[1,1].set_zlabel("density")
    fig.suptitle("Collision rate = " + str(desired_collision_rate))
    fig.tight_layout()
    plt.subplots_adjust(top=0.88)



    # fig, ax = plt.subplots(2,2)
    # c = ax[0,0].imshow(neutral_density_exact_reshaped, extent=[0, Lx, 0, Ly])
    # ax[0,0].set_title("Neutral particle density, pure MC")
    # ax[0,0].set_xlabel("x")
    # ax[0,0].set_ylabel("y")
    # fig.colorbar(c, ax=ax[0,0])

    # c = ax[0,1].imshow(neutral_density_reshaped, extent=[0, Lx, 0, Ly])
    # ax[0,1].set_title("Neutral particle density, KD")
    # ax[0,1].set_xlabel("x")
    # ax[0,1].set_ylabel("y")
    # fig.colorbar(c, ax=ax[0,1])

    # diff = neutral_density_reshaped - neutral_density_exact_reshaped
    # c = ax[1,0].imshow(diff, extent=[0, Lx, 0, Ly])
    # ax[1,0].set_title("Difference")
    # ax[1,0].set_xlabel("x")
    # ax[1,0].set_ylabel("y")
    # fig.colorbar(c, ax=ax[1,0])


    err = np.linalg.norm(diff) / np.linalg.norm(neutral_density_exact)
    print("Error for r = ", desired_collision_rate," is ", err)



    plt.savefig("outputs/mc_r_" + str(desired_collision_rate) + ".png")

plt.show()