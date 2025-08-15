import os
import eval_mut_comb
import json
import pandas as pd
import re 

df_culling = pd.read_csv('/home/s_alya/research_project/Boltz_scores_culled_2fde_preIE.csv')


for file_name in df_culling['name']:
    basename = file_name
    design_pdb = f'/home/s_alya/research_project/Boltz2_results/2fde_results/boltz_results_{basename}/predictions/{basename}/{basename}_model_0_protein.pdb'
 # Extract pose number
    match = re.search(r'pose(\d+)', file_name)
    pose_number = int(match.group(1)) if match else None

    # Remove suffix with variable number
    #pdb_predict_prefix = re.sub(r'_corrected_\d+_model_0_docked\.pdb$', '', file_name)
    #prediction_json = f'/home/s_alya/research_project/chosen_prots_timed/2fde/2fde/pose_{pose_number}/{pdb_predict_prefix}_corrected_predict_prob.csv'
    prediction_json = f'/home/s_alya/research_project/interface_eng_predictions_without_lig/{basename}/{basename}_predicted_mutations.json'
    predict_mut_dict = json.load(open(prediction_json, 'r'))
################################# start from here????????
    
    # 3. get sequence need to engineer
    wt_seq = eval_mut_comb.get_single_chain_sequence(design_pdb)
    mutation_options = eval_mut_comb.filter_mutations(predict_mut_dict,min_ddg=5, max_ddg=20) #####

    pocket_residues = json.load(open(f'/home/s_alya/research_project/chosen_prots_timed/2fde/2fde/pose_{pose_number}/pocket_residues.json', 'r'))

    # use name as identifier e.g. '1flm_HCY4_230_model_0', use 1flm_HCY4 as the prefix
    name_prefix = design_pdb.replace("corrected_", "")
    name_prefix = name_prefix.replace("_model_0_protein.pdb", "")
    name_prefix = name_prefix.replace("boltz_results_", "")
    save_path = f'/home/s_alya/research_project/{name_prefix}_mut_comb.json'

    redesign_residues = [] 
    # add all designed chain, this step is to include all redesigned residues and keep it 
    for chain_id, residues in pocket_residues.items():
        if chain_id.startswith(name_prefix):
            redesign_residues.extend(residues)
    # reindexing the redesign residues to match the mutation options, as chainB or C would add length A
    for i in range(len(redesign_residues)):
        
        if redesign_residues[i] > len(wt_seq):
            # minus the length of single chain
            redesign_residues[i] = int(redesign_residues[i]) - len(wt_seq)
    print(f"Reindexed redesigned residues: {redesign_residues}")

    #4. filter the mutation options we have designed around pockets
    filtered_mutations = [
        mut for mut in mutation_options
        if all(pos not in redesign_residues for pos, aa in mut)
    ]
    print(f"Filtered mutations: {filtered_mutations}")  
       
    all_valid_combinations = eval_mut_comb.generate_combinations(
        mutation_options=mutation_options,
        min_mutations=len(mutation_options)-1,
        max_mutations=len(mutation_options),
        num_combinations=20 #changed to 20
    )
    #generate new fasta file if write_fasta is True 

    eval_mut_comb.generate_fasta(all_valid_combinations,
                                   wt_seq,
                                   save_path,
                                    write_fasta= True)