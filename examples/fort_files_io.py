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

    # These vertices are listed in the order: poloidal, radial, however,
    # the monte carlo code expects them in the order: radial, poloidal. So we need to swap the two columns.
    vertices = np.array(vertices)

    return vertices

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

# Read neighbor list from a fort.35 file
def read_neighbor_list(filename):
    with open(filename, "r") as f:
        # Skip the first line (it is the number of triangles)
        s = f.readline()
        ncells = int(s.split()[0])

        neighbor_list = np.zeros((ncells, 3), dtype=int)
        boundary_conditions = np.zeros((ncells, 3), dtype=int)
        for line in f:
            s = line.split()
            cell_index = int(s[0]) - 1
            neighbor_list[cell_index, :] = [int(s[1]) - 1, int(s[4]) - 1, int(s[7]) - 1]
            boundary_conditions[cell_index, :] = [int(s[3]), int(s[6]), int(s[9])]
            # Not really sure what this next one is, might be surface id in next cell
            surf_id_next_cell = [int(s[2]), int(s[5]), int(s[8])]

    return np.array(neighbor_list), np.array(boundary_conditions)

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

# Read the magnetic field direction from a b2fgmtry file
def read_magnetic_field(filename):
    with open(filename, "r") as f:
        active = False
        magnetic_field = []
        for line in f:
            line = line.strip()
            if line.startswith("*cf:"):
                parts = line.split()
                tag = parts[-1]
                if tag == "cvEb":
                    active = True
                else:
                    active = False
            elif active:
                numbers = list(map(float, line.split()))
                magnetic_field.extend(numbers)

    # TODO: verify the order that cvEb is written in the file
    magnetic_field = np.array(magnetic_field).reshape(3, -1)  # Reshape to (N, 3) for x, y, z components
    magnetic_field = magnetic_field.T  # Transpose to get shape (N, 3)
    for i in range(len(magnetic_field)):
        norm = np.linalg.norm(magnetic_field[i, :])
        if norm != 0:
            magnetic_field[i, :] /= norm  # Normalize the vector
    return magnetic_field


# Using a fort.35 file, convert the plasma background from the
# plasma solver grid to the neutral grid used in the Monte Carlo code
def convert_background_to_neutral_grid(background, magnetic_field, filename):

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
                if species != "e":
                    converted_background[species]["vel_x"][index_neutral_grid] = background[species]["vel"][index_plasma_grid] * magnetic_field[index_plasma_grid][0]
                    converted_background[species]["vel_y"][index_neutral_grid] = background[species]["vel"][index_plasma_grid] * magnetic_field[index_plasma_grid][1]
                    converted_background[species]["vel_z"][index_neutral_grid] = background[species]["vel"][index_plasma_grid] * magnetic_field[index_plasma_grid][2]

    return converted_background

# Read the eirene output from a fort.46 file
# This output is given in the plasma solver grid,
# so we also convert it to the neutral grid using the fort.35 file
def read_eirene_output(filename_46, filename_35):
    # Arrays to hold the data
    atom_density = []
    atom_vel_x = []
    atom_vel_y = []
    atom_vel_z = []
    atom_energy = []

    current_array = None
    with open(filename_46, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("*eirene data field"):
                parts = line.split()
                tag = parts[3]
                
                if tag == "pdena":
                    current_array = atom_density
                elif tag == "vxdena":
                    current_array = atom_vel_x
                elif tag == "vydena":
                    current_array = atom_vel_y
                elif tag == "vzdena":
                    current_array = atom_vel_z
                elif tag == "edena":
                    current_array = atom_energy
                else:
                    current_array = None  # Ignore other types
            
            elif current_array is not None:  # Read numeric values
                numbers = list(map(float, line.split()))
                current_array.extend(numbers)

    # Convert the lists to numpy arrays
    atom_density = np.array(atom_density)
    atom_vel_x = np.array(atom_vel_x)
    atom_vel_y = np.array(atom_vel_y)
    atom_vel_z = np.array(atom_vel_z)
    atom_energy = np.array(atom_energy)

    # Convert from cm^-3 to m^-3
    atom_density *= 1e6
    atom_vel_x *= 1e6
    atom_vel_y *= 1e6
    atom_vel_z *= 1e6
    atom_energy *= 1e6

    # Now we need to convert the eirene output from the plasma grid to the neutral grid
    with open(filename_35, "r") as f:
        # Read the number of points in the neutral grid
        num_points = int(f.readline().strip())
        # Create a dictionary to hold the converted eirene output
        converted_eirene_output = {
            "density": np.zeros(num_points),
            "vel_x": np.zeros(num_points),
            "vel_y": np.zeros(num_points),
            "vel_z": np.zeros(num_points),
            "energy": np.zeros(num_points),
        }

        index_neutral_grid = -1
        for line in f:
            index_neutral_grid+=1
            s = line.split()

            index_plasma_grid = int(s[-4]) - 1

            # If the index is negative, this is a void cell in the neutral grid, so we skip it
            if index_plasma_grid < 0:
                continue

            converted_eirene_output["density"][index_neutral_grid] = atom_density[index_plasma_grid]
            converted_eirene_output["vel_x"][index_neutral_grid] = atom_vel_x[index_plasma_grid]
            converted_eirene_output["vel_y"][index_neutral_grid] = atom_vel_y[index_plasma_grid]
            converted_eirene_output["vel_z"][index_neutral_grid] = atom_vel_z[index_plasma_grid]
            converted_eirene_output["energy"][index_neutral_grid] = atom_energy[index_plasma_grid]

    return converted_eirene_output