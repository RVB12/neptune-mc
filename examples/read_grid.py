import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
import json
import time
import eigrene

vertex_file = "examples/fort.33"
triangle_file = "examples/fort.34"

# Read the vertices
vertices = []
with open(vertex_file, "r") as f:
    # Skip the first line (it is the number of vertices)
    f.readline()
    for line in f:
        # Skip the first number (it is the vertex number)
        s = line.split()
        vertices.append([float(x)/100 for x in s[1:]])


triangles = []
with open(triangle_file, "r") as f:
    # Skip the first line (it is the number of triangles)
    f.readline()
    for line in f:
        # Skip the first number (it is the triangle number)
        s = line.split()
        triangles.append([int(x) - 1 for x in s[1:]])

# print("Vertices:")
# for vertex in vertices:
#     print(vertex)

# print("Triangles:")
# for triangle in triangles:
#     print(triangle)

# Create the triangulation
x = np.array([v[0] for v in vertices])
y = np.array([v[1] for v in vertices])

ncells = len(triangles)
xy = np.concatenate((x[:, None], y[:, None]), axis=1)

config = {}

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)


for i in range(ncells):
    plasma_dens[i] = 1.0e19

for i in range(ncells):
    plasma_temp[i] = 1.0

for i in range(ncells):
    plasma_vel_x[i] = 0

for i in range(ncells):
    plasma_vel_y[i] = 0


config["background"] = {
    "plasma_dens": plasma_dens.tolist(),
    "plasma_temp": plasma_temp.tolist(),
    "plasma_vel_x": plasma_vel_x.tolist(),
    "plasma_vel_y": plasma_vel_y.tolist()
}

config["mesh"] = {
    "type": "triangle",
    "vertices": xy.tolist(),
    "triangles": np.array(triangles).tolist()
}


config["sources"] = [
    {
        "type": "point",
        "n_particles": 1000,
        "intensity": 1.0e10,
        "pos": [0.55, -0.3],
        "vel_mean": [2, 0],
        "vel_stddev": 100,
    },
]

config["boundary_conditions"] = [
    "reflective"
]

config["max_steps"] = 10000000
config["seed"] = 0


inputstring = json.dumps(config)

start_time = time.time()
outputstring = eigrene.simulate(inputstring)
end_time = time.time()

print("Simulation took", end_time - start_time, "seconds")


output = json.loads(outputstring)

stats = output["stats"]
print(stats)

neutral_density = np.array(output["density"])

#
fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor(x, y, neutral_density, triangles=triangles, shading='flat')
fig1.colorbar(tpc)
ax1.set_title('neutral density')

plt.show()

