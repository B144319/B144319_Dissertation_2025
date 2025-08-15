import os
import json
import pandas as pd
import re
import eval_mut_comb
from collections import defaultdict

# Read culled mutation predictions
culling_path = '/home/s_alya/research_project/Boltz_scores_culled_4dri_preIE.csv'
df_culling = pd.read_csv(culling_path)

# Chain lengths for heterodimer (adjust as needed)
chain_lengths = {'A': 121, 'B': 96}

def get_chain_boundaries(chain_lengths):
    chain_boundaries = {}
    start = 1
    for chain_id, length in chain_lengths.items():
        end = start + length - 1
        chain_boundaries[chain_id] = (start, end)
        start = end + 1
    return chain_boundaries

chain_boundaries = get_chain_boundaries(chain_lengths)

def group_mutations_by_position(mutations):
    grouped = defaultdict(list)
    for mut in mutations:
        chain, resnum, aa = mut
        grouped[(chain, resnum)].append(mut)
    return list(grouped.values())

# Process each Boltz design
for file_name in df_culling['name']:
    basename = file_name
    design_pdb = f'/home/s_alya/research_project/Boltz2_results/4dri_results/boltz_results_{basename}/predictions/{basename}/{basename}_model_0_protein.pdb'

    match = re.search(r'pose(\d+)', file_name)
    pose_number = int(match.group(1)) if match else None

    prediction_json = f'/home/s_alya/research_project/interface_eng_predictions_without_lig/{basename}/{basename}_predicted_mutations.json'
    predict_mut_dict = json.load(open(prediction_json, 'r'))

    # Get WT sequences for heterodimer chains
    
    wt_sequences = eval_mut_comb.get_chain_sequences(design_pdb)

    print("Wild-type sequences:", wt_sequences)

    # Filter mutations based on ddG range
    mutation_options = eval_mut_comb.filter_mutations(predict_mut_dict, min_ddg=5, max_ddg=20)
    if any(isinstance(i, list) for i in mutation_options):
        mutation_options = [item for sublist in mutation_options for item in sublist]

    parsed_mutations = []
    print("mutation_options contents:", mutation_options)
    for resnum, aa in mutation_options:
        for chain_id, (start, end) in chain_boundaries.items():
            if start <= resnum <= end:
                resnum_in_chain = resnum - start + 1
                parsed_mutations.append((chain_id, resnum_in_chain, aa))
                break

    pocket_path = f'/home/s_alya/research_project/chosen_prots_timed/4dri/pose_{pose_number}/pocket_residues.json'
    pocket_residues = json.load(open(pocket_path, 'r'))

    name_prefix = basename.replace("corrected_", "").replace("_model_0_protein.pdb", "")
    save_path = f'/home/s_alya/research_project/{name_prefix}'
    os.makedirs(save_path, exist_ok=True)

    redesign_residues = []
    for chain_id, residues in pocket_residues.items():
        redesign_residues.extend([(chain_id[-1], int(res)) for res in residues])

    print(f"Redesigned residues (chain-specific): {redesign_residues}")

    filtered_mutations = [
        mut for mut in parsed_mutations
        if (mut[0], mut[1]) not in redesign_residues
    ]

    print(f"Filtered mutations: {filtered_mutations}")
    grouped_mutations = group_mutations_by_position(filtered_mutations)

    if len(grouped_mutations) == 0:
        print("No mutations passed the filter — using wild-type sequence.")
        all_valid_combinations = [[(chain, pos + 1, aa)
                                   for chain, seq in wt_sequences.items()
                                   for pos, aa in enumerate(seq)]]
    else:
        all_valid_combinations = eval_mut_comb.generate_combinations(
            mutation_options=grouped_mutations,
            min_mutations=len(grouped_mutations) - 1,
            max_mutations=len(grouped_mutations),
            num_combinations=50
        )
        if not all_valid_combinations:
            print("No combinations generated — using wild-type sequence.")
            all_valid_combinations = [[(chain, pos + 1, aa)
                                       for chain, seq in wt_sequences.items()
                                       for pos, aa in enumerate(seq)]]

    # Generate FASTA and comb_info.json
    eval_mut_comb.generate_fasta(
        combinations=all_valid_combinations,
        wt_sequences=wt_sequences,
        save_path=save_path,
        write_fasta=True
    )
