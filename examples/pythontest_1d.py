import numpy as np
import eiron2
from matplotlib import pyplot as plt
import math

Lx = 1
Ly = 1
Lz = 1

nx = 400
ny = 1
nz = 1

ncells = nx * ny * nz

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)
plasma_vel_z = np.zeros(ncells, dtype=np.float64)

for i in range(ncells):
    plasma_dens[i] = 1.0e13

for i in range(ncells):
    plasma_temp[i] = 1.0e5

for i in range(ncells):
    plasma_vel_x[i] = 1.0e2

mesh = eiron2.VoxelMesh(Lx, Ly, Lz, nx, ny, nz)

background = eiron2.PlasmaBackground(
    eiron2.VectorDouble(plasma_dens),
    eiron2.VectorDouble(plasma_temp),
    eiron2.VectorDouble(plasma_vel_x),
    eiron2.VectorDouble(plasma_vel_y),
    eiron2.VectorDouble(plasma_vel_z)
)

params = eiron2.SimulationParams()
params.n_particles = 100000
params.max_steps = 1000000

neutral_density = np.zeros(ncells, dtype=np.float64)

eiron2.simulate(mesh, background, eiron2.VectorDouble(neutral_density), params)

x = np.linspace(0,1,nx)

plt.plot(x, neutral_density)
plt.title("Neutral particle density")
plt.xlabel("x")
plt.show()