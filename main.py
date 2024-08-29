import sys
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import nglview as nv
import matplotlib.pyplot as plt

def get_input():
    if len(sys.argv) != 3:
        print("Usage: python script.py <protein_pdb> <ligand_pdb>")
        sys.exit(1)

    protein_pdb = sys.argv[1]
    ligand_pdb = sys.argv[2]

    protein_u = mda.Universe(protein_pdb)
    ligand_u = mda.Universe(ligand_pdb)

    return protein_u, ligand_u

def process_data(protein_u, ligand_u):
    protein = protein_u.select_atoms("protein")
    ligand = ligand_u.select_atoms("all")

    ligand.translate(protein.center_of_mass() - ligand.center_of_mass())

    distances_result = distances.distance_array(protein.positions, ligand.positions)
    min_distances = distances_result.min(axis=1)

    return protein_u, ligand_u, min_distances

def display_output(protein_u, ligand_u, min_distances):
    view = nv.show_mdanalysis(protein_u)
    view.add_representation("cartoon", selection="protein")
    view.add_component(ligand_u)
    view.add_representation("ball+stick", component=1)
    view

    plt.hist(min_distances, bins=50)
    plt.xlabel("Distance (Å)")
    plt.ylabel("Frequency")
    plt.title("Distance Distribution between Protein and Ligand")
    plt.show()

    print(f"Minimum distance: {min_distances.min():.2f} Å")
    print(f"Average distance: {min_distances.mean():.2f} Å")

def main():
    protein_u, ligand_u = get_input()
    protein_u, ligand_u, min_distances = process_data(protein_u, ligand_u)
    display_output(protein_u, ligand_u, min_distances)

if __name__ == "__main__":
    main()