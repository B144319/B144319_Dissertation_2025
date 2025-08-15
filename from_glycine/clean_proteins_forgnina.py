import pymol
from pymol import cmd
import os 

#input and output directory: 
input_dir = 'top_protein_pdbs'
output_dir = 'cleaned_protein_pdbs'
if not os.path.exists(output_dir): 
    os.makedirs(output_dir)

#Initialize pymol: 
pymol.finish_launching()

#Clean each pdb file - remove non-protein atoms e.g. solvents and ligands 
for pdb_file in os.listdir(input_dir): 
    cmd.reinitialize()
    cmd.load(os.path.join(input_dir, pdb_file))
    cmd.remove("hetatm")
    cmd.remove("solvent")
    output_path = os.path.join(output_dir, f"cleaned_{pdb_file}")
    cmd.save(output_path)

    