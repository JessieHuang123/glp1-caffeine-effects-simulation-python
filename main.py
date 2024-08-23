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

# Visualization: Display molecular structures using NGLView
view = nv.show_mdanalysis(protein_u)
view.add_representation("cartoon", selection="protein")
view.add_component(ligand_u)  # Add caffeine
view.add_representation("ball+stick", component=1)  # Display caffeine in ball and stick style
view

# Calculate and visualize minimum distances
min_distances = distances_result.min(axis=1)
plt.hist(min_distances, bins=50)
plt.xlabel("Distance (Å)")
plt.ylabel("Frequency")
plt.title("Distance Distribution between Protein and Caffeine")
plt.show()

# Print results
print(f"Minimum distance: {min_distances.min():.2f} Å")
print(f"Average distance: {min_distances.mean():.2f} Å")