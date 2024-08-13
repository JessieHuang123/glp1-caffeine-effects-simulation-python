import MDAnalysis as mda
from MDAnalysis.analysis import distances
import nglview as nv
import matplotlib.pyplot as plt

sample_pdb = "sample.pdb"   # sample pdb file

# Loading MDAnalysis
protein_u = mda.Universe(sample_pdb)

