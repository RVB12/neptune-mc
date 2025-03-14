import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math
import json
import time
import eigrene

vertex_file = "examples/fort.33"
triangle_file = "examples/fort.34"
mapping_file = "examples/fort.35"
background_file = "examples/b2fstati"
geometry_file = "examples/b2fgmtry"

# Read the vertices
vertices = []
with open(vertex_file, "r") as f:
    # Skip the first line (it is the number of vertices)
    f.readline()
    for line in f:
        # Skip the first number (it is the vertex number)
        s = line.split()
        vertices.append([float(x)/100 for x in s[1:]])

# Read the triangles
triangles = []
with open(triangle_file, "r") as f:
    # Skip the first line (it is the number of triangles)
    f.readline()
    for line in f:
        # Skip the first number (it is the triangle number)
        s = line.split()
        triangles.append([int(x) - 1 for x in s[1:]])

# Read the magnetic field from the b2fgmtry file
magnetic_field_plasma_grid = []
current_array = None
with open(geometry_file, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith("*cf:"):
            parts = line.split()
            tag = parts[-1]
            
            if tag == "cvEb":
                current_array = magnetic_field_plasma_grid
            else:
                current_array = None  # Ignore other types
        
        elif current_array is not None:  # Read numeric values
            numbers = list(map(float, line.split()))
            current_array.extend(numbers)

# The magnetic field array is actually the magnetic field in the x, y, and z directions
# reshape it to be a list of 3-tuples
magnetic_field_plasma_grid = np.array(magnetic_field_plasma_grid).reshape(-1, 3)


# Read the densities from the b2fstati file
ion_density_plasma_grid = []
electron_density_plasma_grid = []
ion_temp_plasma_grid = []
ion_vel_parallel_plasma_grid = []
current_array = None
with open(background_file, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith("*cf:"):
            parts = line.split()
            tag = parts[-1]
            
            if tag == "na":
                current_array = ion_density_plasma_grid
            elif tag == "ne":
                current_array = electron_density_plasma_grid
            elif tag == "ti":
                current_array = ion_temp_plasma_grid
            elif tag == "ua":
                current_array = ion_vel_parallel_plasma_grid
            else:
                current_array = None  # Ignore other types
        
        elif current_array is not None:  # Read numeric values
            numbers = list(map(float, line.split()))
            current_array.extend(numbers)

# The density array is actually the density of the atoms and then the density of the ions
# only keep the latter half
ion_density_plasma_grid = ion_density_plasma_grid[len(ion_density_plasma_grid)//2:]

print("ion density")
print(ion_density_plasma_grid[0:10])

# ion temp is defined in joules but we need it in eV
ion_temp_plasma_grid = [x / 1.602176634e-19 for x in ion_temp_plasma_grid]

print("ion temp")
print(ion_temp_plasma_grid[0:10])
print(max(ion_temp_plasma_grid))

# Create the triangulation
x = np.array([v[0] for v in vertices])
y = np.array([v[1] for v in vertices])

ncells = len(triangles)
xy = np.concatenate((x[:, None], y[:, None]), axis=1)

plasma_dens = np.zeros(ncells, dtype=np.float64)
plasma_temp = np.zeros(ncells, dtype=np.float64)
plasma_vel_x = np.zeros(ncells, dtype=np.float64)
plasma_vel_y = np.zeros(ncells, dtype=np.float64)
plasma_vel_z = np.zeros(ncells, dtype=np.float64)
plasma_vel_parallel = np.zeros(ncells, dtype=np.float64)

# Read the mapping
ind = -1
with open(mapping_file, "r") as f:
    # Skip the first line (it is the number of triangles)
    f.readline()
    for line in f:
        ind+=1
        s = line.split()
        plasma_cell_index = int(s[-4]) - 1

        if(plasma_cell_index < 0):
            continue

        magnetic_unit = magnetic_field_plasma_grid[plasma_cell_index,:]
        magnetic_unit = magnetic_unit / np.linalg.norm(magnetic_unit)


        plasma_dens[ind] = ion_density_plasma_grid[plasma_cell_index]
        plasma_temp[ind] = ion_temp_plasma_grid[plasma_cell_index]
        plasma_vel_parallel[ind] = ion_vel_parallel_plasma_grid[plasma_cell_index]

        plasma_vel_x[ind] = ion_vel_parallel_plasma_grid[plasma_cell_index] * magnetic_unit[0]
        plasma_vel_y[ind] = ion_vel_parallel_plasma_grid[plasma_cell_index] * magnetic_unit[1]
        plasma_vel_z[ind] = ion_vel_parallel_plasma_grid[plasma_cell_index] * magnetic_unit[2]




# print("Vertices:")
# for vertex in vertices:
#     print(vertex)

# print("Triangles:")
# for triangle in triangles:
#     print(triangle)


config = {}



# for i in range(ncells):
#     plasma_dens[i] = 1.0e19

# for i in range(ncells):
#     plasma_temp[i] = 1.0

# for i in range(ncells):
#     plasma_vel_x[i] = 0

# for i in range(ncells):
#     plasma_vel_y[i] = 0


config["background"] = {
    "plasma_dens": plasma_dens.tolist(),
    "plasma_temp": plasma_temp.tolist(),
    "plasma_vel_x": plasma_vel_x.tolist(),
    "plasma_vel_y": plasma_vel_y.tolist(),
    "plasma_vel_z": plasma_vel_z.tolist(),
}

config["mesh"] = {
    "type": "triangle",
    "vertices": xy.tolist(),
    "triangles": np.array(triangles).tolist()
}

config["sources"] = [
    {
        "type": "point",
        "n_particles": 100000,
        "intensity": 1.0e19,
        "pos": [0.441, 0.0, -0.2],
        "vel_mean": [0, 0., 0],
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
tpc = ax1.tripcolor(x, y, np.log10(plasma_dens+1), triangles=triangles, shading='flat', vmin=19, vmax=21)
fig1.colorbar(tpc)
ax1.set_title('ion density (m^-3)')

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor(x, y, plasma_temp, triangles=triangles, shading='flat')
fig1.colorbar(tpc)
ax1.set_title('ion temp (eV)')

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor(x, y, np.log10(np.abs(plasma_vel_x) + 1), triangles=triangles, shading='flat')
fig1.colorbar(tpc)
ax1.set_title('background velocity (x)')

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor(x, y, np.log10(np.abs(plasma_vel_y) + 1), triangles=triangles, shading='flat')
fig1.colorbar(tpc)
ax1.set_title('background velocity (y)')

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor(x, y, np.log10(np.abs(plasma_vel_z) + 1), triangles=triangles, shading='flat')
fig1.colorbar(tpc)
ax1.set_title('background velocity (z)')

fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tpc = ax1.tripcolor(x, y, neutral_density, triangles=triangles, shading='flat')
fig1.colorbar(tpc)
ax1.set_title('neutral density (m^-3)')

plt.show()

