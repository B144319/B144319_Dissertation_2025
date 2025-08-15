import os
import csv
import pandas as pd

# Paths
cull_file = "Boltz_scores_culled_2fde_postIE.csv"
input_dir = "Boltz_after_int_eng/2fde_results"
output_dir = "cleaned_prots_2fde_postIE"
os.makedirs(output_dir, exist_ok=True)

# Reading protein names from CSV
df = pd.read_csv(cull_file)

for prot in df['name']:  
    basename = prot
    pdb_path = os.path.join(
        input_dir,
        f"boltz_results_{basename}",
        "predictions",
        basename,
        f"{basename}_model_0.pdb"
    )

    if not os.path.exists(pdb_path):
        continue

    # Output paths
    cleaned_path = os.path.join(output_dir, f"{basename}_cleaned.pdb")
    hetatm_path = os.path.join(output_dir, f"{basename}_ligand.pdb")

    with open(pdb_path, "r") as infile, \
         open(cleaned_path, "w") as clean_out, \
         open(hetatm_path, "w") as hetatm_out:

        for line in infile:
            if line.startswith("HETATM"):
                hetatm_out.write(line)
            else:
                clean_out.write(line)

