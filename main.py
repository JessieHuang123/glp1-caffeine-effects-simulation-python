import sys
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import nglview as nv
import matplotlib.pyplot as plt

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script.py <protein_pdb> <ligand_pdb>")
    sys.exit(1)

# Get the PDB file paths from command-line arguments
protein_pdb = sys.argv[1]
ligand_pdb = sys.argv[2]

# Load structures using MDAnalysis
protein_u = mda.Universe(protein_pdb)
ligand_u = mda.Universe(ligand_pdb)

# Select protein and ligand
protein = protein_u.select_atoms("protein")
ligand = ligand_u.select_atoms("all")  # Select all atoms of ligand

# Move ligand to the center of the protein (manual adjustment may be needed in this case)
ligand.translate(protein.center_of_mass() - ligand.center_of_mass())

# Calculate distances between protein and ligand
distances_result = distances.distance_array(protein.positions, ligand.positions)

# Visualization: Display molecular structures using NGLView
view = nv.show_mdanalysis(protein_u)
view.add_representation("cartoon", selection="protein")
view.add_component(ligand_u)  # Add ligand
view.add_representation("ball+stick", component=1)  # Display ligand in ball and stick style
view

# Calculate and visualize minimum distances
min_distances = distances_result.min(axis=1)
plt.hist(min_distances, bins=50)
plt.xlabel("Distance (Å)")
plt.ylabel("Frequency")
plt.title("Distance Distribution between Protein and Ligand")
plt.show()

# Print results
print(f"Minimum distance: {min_distances.min():.2f} Å")
print(f"Average distance: {min_distances.mean():.2f} Å")
