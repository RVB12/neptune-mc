import h5py
import numpy as np

f = h5py.File("build/input.h5", "w")

n_cells = 10

# Background
background = f.create_group("background")
Dplus = background.create_group("D+")
Dplus["density"] = 1.0e21 * np.ones(n_cells, dtype=np.float64)
Dplus["density"].attrs["units"] = "m^-3"
Dplus["velocity_x"] = np.zeros(n_cells, dtype=np.float64)
Dplus["velocity_x"].attrs["units"] = "u*m/s"
Dplus["velocity_y"] = np.zeros(n_cells, dtype=np.float64)
Dplus["velocity_y"].attrs["units"] = "u*m/s"
Dplus["velocity_z"] = np.zeros(n_cells, dtype=np.float64)
Dplus["velocity_z"].attrs["units"] = "u*m/s"
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

sources = f.create_group("sources")

sources["n_particles"] = 10000
sources["stratification"] = "source_strength"

# Some volume recombination sources
source = sources.create_group("R_vol")
source["type"] = "recombination"
source["background_species"] = "D+"
source["particle_species"] = "D"
source["cell_indices"] = np.arange(n_cells, dtype=np.int64)

# Also some wall recycling sources
source = sources.create_group("R_wall")
source["type"] = "wall_recycling"
source["background_species"] = "D+"
source["particle_species"] = "D"
source["cell_indices"] = np.array([n_cells - 1], dtype=np.int64)
source["surface_indices"] = np.array([3], dtype=np.int64)
source["sheath_potential"] = 3.1
source["franck_condon_temp"] = 2.0
source["franck_condon_temp"].attrs["units"] = "eV"
source["Rt"] = 0.9
source["Rr"] = 0.9

# Estimators
estimators = f.create_group("estimators")
D_estimator = estimators.create_group("D")
D_estimator["species"] = "D"
D_estimator["qoi"] = "density"
D_estimator["estimator_type"] = "track_length"
D_estimator["estimate_variance"] = 0

# Boundary conditions
boundary_conditions = f.create_group("boundary_conditions")

for i in range(6):
    bc = boundary_conditions.create_group(f"bc{i}")
    bc["type"] = "specular"

