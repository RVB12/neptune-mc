import numpy as np

# Read vertices from a fort.33 file
def read_vertices(filename):
    vertices = []
    with open(filename, "r") as f:
        # Skip the first line (it is the number of vertices)
        f.readline()
        for line in f:
            # Skip the first number (it is the vertex number)
            s = line.split()
            vertices.append([float(x)/100 for x in s[1:]])

    return np.array(vertices)

# Read triangles from a fort.34 file
def read_triangles(filename):
    triangles = []
    with open(filename, "r") as f:
        # Skip the first line (it is the number of triangles)
        f.readline()
        for line in f:
            # Skip the first number (it is the triangle number)
            s = line.split()
            triangles.append([int(x) - 1 for x in s[1:]])

    return np.array(triangles)

# Convert species information from a list of atomic mass and charge to species names
def convert_species(zamin, zamax, nuclear_charge, atomic_mass):
    species = []
    for i in range(len(zamin)):
        if zamin[i] == 0.0 and zamax[i] == 0.0 and nuclear_charge[i] == 1.0 and atomic_mass[i] == 1.0:
            species.append("H")
        elif zamin[i] == 1.0 and zamax[i] == 1.0 and nuclear_charge[i] == 1.0 and atomic_mass[i] == 1.0:
            species.append("H+")
        elif zamin[i] == 0.0 and zamax[i] == 0.0 and nuclear_charge[i] == 1.0 and atomic_mass[i] == 2.0:
            species.append("D")
        elif zamin[i] == 1.0 and zamax[i] == 1.0 and nuclear_charge[i] == 1.0 and atomic_mass[i] == 2.0:
            species.append("D+")
    return species

# Read information about the plasma background from a b2fstati file
# This is given in the grid of the plasma solver
def read_background(filename):

    # Arrays to hold the data
    atom_density = []
    atom_vel = []
    neutral_temp = []
    ion_temp = []
    electron_density = []
    electron_temp = []
    # electron_vel = [] # not used in the code

    # Arrays to hold species info
    zamin = []
    zamax = []
    nuclear_charge = []
    atomic_mass = []

    current_array = None
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("*cf:"):
                parts = line.split()
                tag = parts[-1]
                
                if tag == "zamin":
                    current_array = zamin
                elif tag == "zamax":
                    current_array = zamax
                elif tag == "zn":
                    current_array = nuclear_charge
                elif tag == "am":
                    current_array = atomic_mass
                elif tag == "na":
                    current_array = atom_density
                elif tag == "ne":
                    current_array = electron_density
                elif tag == "ti":
                    current_array = ion_temp
                elif tag == "tn":
                    current_array = neutral_temp
                elif tag == "te":
                    current_array = electron_temp
                elif tag == "ua":
                    current_array = atom_vel
                # elif tag == "fna":
                #     current_array = ion_flux_plasma_grid
                else:
                    current_array = None  # Ignore other types
            
            elif current_array is not None:  # Read numeric values
                numbers = list(map(float, line.split()))
                current_array.extend(numbers)

    # Convert species information to species names
    species = convert_species(zamin, zamax, nuclear_charge, atomic_mass)

    atom_density = np.array_split(np.array(atom_density), len(species))
    atom_vel = np.array_split(np.array(atom_vel), len(species))
    neutral_temp = np.array(neutral_temp)
    ion_temp = np.array(ion_temp)

    background = {}
    for i in range(len(species)):
        background[species[i]] = {
            "density": atom_density[i],
            "vel": atom_vel[i],
        }
        if species[i] == "H" or species[i] == "D":
            background[species[i]]["temp"] = neutral_temp
        else:
            background[species[i]]["temp"] = ion_temp

    background["e"] = {
        "density": np.array(electron_density),
        "temp": np.array(electron_temp),
    }

    # Convert the temperature from joules to eV
    for species in background.keys():
        background[species]["temp"] = background[species]["temp"] / 1.602176634e-19

    return background

# Using a fort.35 file, convert the plasma background from the
# plasma solver grid to the neutral grid used in the Monte Carlo code
def convert_background_to_neutral_grid(background, filename):

    with open(filename, "r") as f:
        # Read the number of points in the neutral grid
        num_points = int(f.readline().strip())
        # Create a dictionary to hold the converted background
        converted_background = {}
        for species in background.keys():
            converted_background[species] = {
                "density": np.zeros(num_points),
                "temp": np.zeros(num_points),
            }
            if species != "e":
                converted_background[species]["vel_x"] = np.zeros(num_points)
                converted_background[species]["vel_y"] = np.zeros(num_points)
                converted_background[species]["vel_z"] = np.zeros(num_points)


        index_neutral_grid = -1
        for line in f:
            index_neutral_grid+=1
            s = line.split()

            index_plasma_grid = int(s[-4]) - 1

            # If the index is negative, this is a void cell in the neutral grid, so we skip it
            if index_plasma_grid < 0:
                continue

            for species in background.keys():
                converted_background[species]["density"][index_neutral_grid] = background[species]["density"][index_plasma_grid]
                converted_background[species]["temp"][index_neutral_grid] = background[species]["temp"][index_plasma_grid]
                # TODO: Convert parallel velocity to x, y, z components.

    return converted_background
