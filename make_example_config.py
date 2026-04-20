import h5py
import numpy as np

f = h5py.File("input.h5", "w")

n_cells = 10

# Background
background = f.create_group("background")
Dplus = background.create_group("Dplus")
Dplus["density"] = np.ones(n_cells, dtype=np.float64)
Dplus["density"].attrs["units"] = "m^-3"
Dplus["momentum_x"] = np.zeros(n_cells, dtype=np.float64)
Dplus["momentum_x"].attrs["units"] = "u*m/s"
Dplus["momentum_y"] = np.zeros(n_cells, dtype=np.float64)
Dplus["momentum_y"].attrs["units"] = "u*m/s"
Dplus["momentum_z"] = np.zeros(n_cells, dtype=np.float64)
Dplus["momentum_z"].attrs["units"] = "u*m/s"
Dplus["temperature"] = np.ones(n_cells, dtype=np.float64)
Dplus["temperature"].attrs["units"] = "eV"

# Mesh
mesh = f.create_group("mesh")

mesh["type"] = "uniform"
mesh["Lx"] = 1.0
mesh["Lx"].attrs["units"] = "m"
mesh["Ly"] = 1.0
mesh["Ly"].attrs["units"] = "m"
mesh["Lz"] = 1.0
mesh["Lz"].attrs["units"] = "m"
mesh["nx"] = n_cells
mesh["ny"] = n_cells
mesh["nz"] = n_cells

# Collisions

collisions = f.create_group("collisions")

coll0 = collisions.create_group("0")
coll0["background_species"] = "D+"
coll0["particle_species"] = "D"
coll0["type"] = "charge_exchange"
coll0["model"] = "simple"

coll1 = collisions.create_group("1")
coll1["background_species"] = "D+"
coll1["particle_species"] = "D"
coll1["type"] = "ionisation"
coll1["model"] = "simple"

# Sources




