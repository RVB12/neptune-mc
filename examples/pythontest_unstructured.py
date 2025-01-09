import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import eiron2

# First create the x and y coordinates of the points.
# n_angles = 36
# n_radii = 8
# min_radius = 0.25
# radii = np.linspace(min_radius, 0.95, n_radii)

# angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
# angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
# angles[:, 1::2] += np.pi / n_angles

# x = (radii * np.cos(angles)).flatten()
# y = (radii * np.sin(angles)).flatten()
# z = (np.cos(radii) * np.cos(3 * angles)).flatten()

# x = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 10))[0].flatten()
# y = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 10))[1].flatten()

x = np.array([0., 2., 0.])
y = np.array([0., 0., 2.])

# Create the Triangulation; no triangles so Delaunay triangulation created.
triang = tri.Triangulation(x, y)

# Mask off unwanted triangles.
# triang.set_mask(np.hypot(x[triang.triangles].mean(axis=1),
#                          y[triang.triangles].mean(axis=1))
#                 < min_radius)


ncells = triang.triangles.shape[0]
xy = np.concatenate((x[:, None], y[:, None]), axis=1).flatten(order='C')
triangles = triang.triangles.flatten(order='C')

mesh = eiron2.TriangleMesh(eiron2.VectorDouble(xy), eiron2.VectorInt(triangles))

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)

for i in range(ncells):
    plasma_dens[i] = 1.0e13
    plasma_vel_x[i] = 10
    plasma_vel_y[i] = 2
    plasma_temp[i] = 0.25

background = eiron2.PlasmaBackground(
    eiron2.VectorDouble(plasma_dens),
    eiron2.VectorDouble(plasma_temp),
    eiron2.VectorDouble(plasma_vel_x),
    eiron2.VectorDouble(plasma_vel_y)
)

params = eiron2.SimulationParams()
params.n_particles = 1
params.max_steps = 1000

neutral_density = np.zeros(ncells, dtype=np.float64)

eiron2.simulate(mesh, background, eiron2.VectorDouble(neutral_density), params)


fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor(triang, neutral_density, shading='flat')
fig1.colorbar(tpc)
ax1.set_title('neutral density')

plt.show()