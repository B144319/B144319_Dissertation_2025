import eval_mut_comb
import json
import pandas as pd 
import re 

# loop through the culling df.prot_name:

df_culling = pd.read_csv('/home/s_alya/research_project/gnina_cull_2fde.csv')
import os
import re
import pandas as pd
import eval_mut_comb

df_culling = pd.read_csv('/home/s_alya/research_project/gnina_cull_2fde.csv')

for file_name in df_culling['filename']:
    basename = file_name.replace("_model_0_docked.pdb", "")
    pdb_path = f'/home/s_alya/research_project/Boltz2_results/2fde_results/boltz_results_{basename}/predictions/{basename}/{basename}_combined.pdb'

    # Extract pose number
    match = re.search(r'pose(\d+)', file_name)
    pose_number = int(match.group(1)) if match else None

    # Remove suffix with variable number
    pdb_predict_prefix = re.sub(r'_corrected_\d+_model_0_docked\.pdb$', '', file_name)

    redesign_path = f'/home/s_alya/research_project/chosen_prots_timed/2fde/2fde/pose_{pose_number}/{pdb_predict_prefix}_corrected_predict_prob.csv'
    output_prefix = basename
    output_json = f"{output_prefix}/{output_prefix}_predicted_mutations.json"

    # Skip if output already exists
    if os.path.exists(output_json):
        print(f"Skipping {output_prefix}: already processed.")
        continue

    print(f"Running pipeline for {output_prefix}")
    print(f"Using redesign csv: {redesign_path}")
    eval_mut_comb.do_mut_pipeline(pdb_path, redesign_path, output_prefix)
