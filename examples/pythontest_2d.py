import numpy as np
import eiron2
from matplotlib import pyplot as plt
import math

Lx = 1
Ly = 1

nx = 400
ny = 400

ncells = nx * ny

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)

for i in range(nx):
    for j in range(ny):
        # plasma_dens[i + j*nx] = 100*(1.0+np.sin(j/100*2*math.pi))
        plasma_dens[i + j*nx] = 1.0e13

for i in range(ncells):
    plasma_vel_x[i] = 10
    plasma_vel_y[i] = 2

for i in range(ncells):
    plasma_temp[i] = 0.25

mesh = eiron2.VoxelMesh(Lx, Ly, nx, ny)

background = eiron2.PlasmaBackground(
    eiron2.VectorDouble(plasma_dens),
    eiron2.VectorDouble(plasma_temp),
    eiron2.VectorDouble(plasma_vel_x),
    eiron2.VectorDouble(plasma_vel_y)
)

params = eiron2.SimulationParams()
params.n_particles = 100000
params.max_steps = 10000000

neutral_density = np.zeros(ncells, dtype=np.float64)

eiron2.simulate(mesh, background, eiron2.VectorDouble(neutral_density), params)

neutral_density_reshaped = np.reshape(neutral_density, (nx, ny))

c = plt.imshow(neutral_density_reshaped, extent=[0, 1, 0, 1])
plt.title("Neutral particle density (log)")
plt.colorbar(c)
plt.xlabel("x")
plt.ylabel("y")
plt.show()