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
