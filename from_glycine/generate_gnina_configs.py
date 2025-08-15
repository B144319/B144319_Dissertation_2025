import os

# Directories
og_lig_dir = "original_ligands"
protein_dir = "cleaned_proteins_pdbs"
myc_rotamer_dir = "mycolactone_rotamers"
config_dir = "gnina_configs"

# Protein names
protein_names = ['1b6k', '1yt9', '2fdd', '2fde', '2i4d', '3o9b', '6d0d', '6dp0e', '6oxp', '6uwc', '6w6s', '6w6t','6dh1', '6dh4', '6dh7', '3o99', '6dgy', '5gpg',
           '1fap', '4drh', '4dri', '4drj', '3o9a', '4m8x',
           '4m8y', '2i4w', '2i4x', '3fap', '4fap', '3mxe',
           '3gi6', '3o9f', '3o9g', '6oxq', '60de']

# List of all rotamer files in rotamer directory
myc_rotamers = sorted([f for f in os.listdir(myc_rotamer_dir) if f.endswith(".sdf")])

# Make config output directory
os.makedirs(config_dir, exist_ok=True)

# Generate config files
for protein in protein_names:
    receptor = os.path.join(protein_dir, f"{protein}_cleaned.pdb")
    autobox_lig = os.path.join(og_lig_dir, f"og_lig_{protein}.mol2")

    for rotamer in myc_rotamers:
        rotamer_name = os.path.splitext(rotamer)[0]
        ligand_path = os.path.join(myc_rotamer_dir, rotamer)
        out_path = f"/home/s_alya/research_project/{protein}_{rotamer_name}_docked.pdb"

        config_content = f"""receptor = /home/s_alya/research_project/{receptor}
ligand = /home/s_alya/research_project/{ligand_path}
autobox_ligand = /home/s_alya/research_project/{autobox_lig}
out = {out_path}
num_modes = 5
exhaustiveness = 8
"""

        config_filename = f"{protein}_{rotamer_name}.txt"
        with open(os.path.join(config_dir, config_filename), "w") as f:
            f.write(config_content)