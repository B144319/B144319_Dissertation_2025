import os
import re
import pandas as pd
from multiprocessing import Pool, cpu_count
import eval_mut_comb

def process_one_protein(file_name):
    """
    Process a single PDB entry by running the mutation pipeline.

    Args:
        file_name (str): Filename from the culling CSV.

    Returns:
        str: Status message
    """
    try:
        basename = file_name.replace("_model_0_docked.pdb", "")
        pdb_path = f'/home/s_alya/research_project/Boltz2_results/4dri_results/boltz_results_{basename}/predictions/{basename}/{basename}_combined.pdb'

        # Extract pose number
        match = re.search(r'pose(\d+)', file_name)
        pose_number = int(match.group(1)) if match else None

        # Remove suffix with variable number
        pdb_predict_prefix = re.sub(r'_corrected_\d+_model_0_docked\.pdb$', '', file_name)

        redesign_path = f'/home/s_alya/research_project/chosen_prots_timed/4dri/pose_{pose_number}/{pdb_predict_prefix}_corrected_predict_prob.csv'
        output_prefix = basename
        output_json = f"{output_prefix}/{output_prefix}_predicted_mutations.json"

        # Skip if already done
        if os.path.exists(output_json):
            return f"Skipping {output_prefix}: already processed."

        print(f"Running pipeline for {output_prefix}")
        print(f"Using redesign csv: {redesign_path}")
        eval_mut_comb.do_mut_pipeline(pdb_path, redesign_path, output_prefix)
        return f"Done: {output_prefix}"

    except Exception as e:
        return f"Error processing {file_name}: {str(e)}"

if __name__ == "__main__":
    df_culling = pd.read_csv('/home/s_alya/research_project/gnina_cull_4dri.csv')
    filenames = df_culling['filename'].tolist()

    # Use 75% of CPU cores
    num_processes = max(1, int(cpu_count() * 0.75))
    print(f"Using {num_processes} CPU cores.")

    with Pool(processes=num_processes) as pool:
        results = pool.map(process_one_protein, filenames)

    # Print the results
    for result in results:
        print(result)

