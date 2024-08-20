import MDAnalysis as mda
from MDAnalysis.analysis import distances
import nglview as nv
import matplotlib.pyplot as plt

# PDB file paths
protein_pdb = "7ki0.pdb"   # Protein PDB file
ligand_pdb = "Caffeine.pdb" # Caffeine PDB file

# Load structures using MDAnalysis
protein_u = mda.Universe(protein_pdb)
ligand_u = mda.Universe(ligand_pdb)

# Select protein and ligand
protein = protein_u.select_atoms("protein")
ligand = ligand_u.select_atoms("all")  # Select all atoms of caffeine

# Move caffeine to the center of the protein (manual adjustment may be needed in this case)
ligand.translate(protein.center_of_mass() - ligand.center_of_mass())

# Calculate distances between protein and caffeine
distances_result = distances.distance_array(protein.positions, ligand.positions)