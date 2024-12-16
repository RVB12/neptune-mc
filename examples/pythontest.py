import numpy as np
import eiron2
from matplotlib import pyplot as plt

Lx = 1
Ly = 1
Lz = 1

nx = 100
ny = 1
nz = 1

ncells = nx * ny * nz

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)
plasma_vel_z = np.zeros(ncells, dtype=np.float64)

for i in range(ncells):
    plasma_dens[i] = 1

for i in range(ncells):
    plasma_temp[i] = 0.25

mesh = eiron2.VoxelMesh(Lx, Ly, Lz, nx, ny, nz)

background = eiron2.PlasmaBackground(
    eiron2.VectorDouble(plasma_dens),
    eiron2.VectorDouble(plasma_temp),
    eiron2.VectorDouble(plasma_vel_x),
    eiron2.VectorDouble(plasma_vel_y),
    eiron2.VectorDouble(plasma_vel_z)
)

params = eiron2.SimulationParams()

neutral_density = np.zeros(ncells, dtype=np.float64)

eiron2.simulate(mesh, background, eiron2.VectorDouble(neutral_density), params)

plt.plot(neutral_density)
plt.show()
