import cotimedutils
import os
import json
import ampal

#Based on cotimedutils.py by Handing Wang. 
#Edited by student B144319

def do_seq_design_pipeline(protein_name: str, pocket_size: float, cutoff: float) -> list:

    '''
        The pipeline for sequence design, including:
        1. Process the generated CSV file from the timed design CSV, 'cnn_timed_50percentmodel_cofactor-148-1.73.csv'
        and 'datasetmap.txt', to get the predicted probabilities for each residue.
        2. Process ligand surrounding residues, 'pocket_residues.json',
        3. Filtering the residues based on the predicted probabilities and the cutoff value. 
        get a JSON file with the predicted sequence for each residue.e.g."1flm_hcy12_redesign_pocket_sequence.json"ArithmeticError
        4. Generate the mutation combinations based on the predicted JSON file.
        and here we assume the protein is a homodimer, so we need to mutate all residues in the single chain.
        Arguments:
            protein_name (str): The name of the protein, e.g. '1flm_hcy12'.
        
        Returns: generates a JSON file with the predicted sequence for each residue.


    '''

    # --- Section 1: Process generated CSV file from timed design CSV ---

    pdb_name = protein_name # e.g. '3cyl_TES_A_124' in the form of pdbid_lig_ligchain_ligid
    probabilities_path = 'CoTIMED.csv'
    label_map_path = 'datasetmap.txt'
    output_csv_name = f'{pdb_name}_predict_prob.csv'

    # Process data and get the DataFrame
    df_prob = cotimedutils.process_data(probabilities_path, label_map_path, output_csv_name)
    
    # --- Section 2: Process the redesign sequence and JSON file from the selected sequence ---
    cotimedutils.get_pocket_residues(pocket_size)
    current_directory = os.getcwd()

    # Load redesign sequence information from a JSON file

    json_input_path = os.path.join(current_directory, 'pocket_residues.json')
    with open(json_input_path, 'r') as f:
        redesign_sequence = json.load(f)

    predict_seq = {}
    for key_full_name in redesign_sequence.keys():

        predict_seq[key_full_name] = {}  # Initialize a sub-dictionary for each key_full_name
        for item_idx in redesign_sequence[key_full_name]:
            # Ensure item_idx is an integer for lookup
            predict_seq[key_full_name][item_idx] = cotimedutils.get_amino_acid(key_full_name, int(item_idx), df_prob, cutoff)

    # Save the predicted sequence to a JSON file
    output_json_name = f'{pdb_name}_redesign_pocket_sequence.json'
    json_output_path = os.path.join(current_directory, output_json_name)

    with open(json_output_path, 'w') as f:
        json.dump(predict_seq, f, indent=4)
    print(f"Predicted sequence data saved to '{output_json_name}'.")

    # --- Section 3: Generate the mutatation combinations based on the predicted JSON file ---
    # Generate the mutation combinations
    process_pdb = protein_name + '.pdb'
    redesign_sequence = cotimedutils.generate_fasta(process_pdb)
    mutation_combinations = cotimedutils.generate_combination(prot_name, json_output_path)
    print(f"Generated {len(mutation_combinations)} mutation combinations.")
    print(f"Mutation combinations: {mutation_combinations}")
    
    # Save the mutation combinations to a CSV file
    #import pandas as pd
    #mutation_combinations = pd.DataFrame(mutation_combinations)
    #mutation_combinations.to_csv(f'{prot_name}_mutation_combinations.csv', index=False)
    # as a heterodimer, both chains need to be mutated. Combine both chains into a single line, mutate then split. 
    
    double_chain = cotimedutils.get_double_chain(process_pdb)
    print(f"Combined chain for {prot_name}: {double_chain}")

    mutated_dict = cotimedutils.all_on_1chain(mutation_combinations, double_chain, prot_name) 
    
    # Split into chain A and chain B in a new dict, using ampal to get the chain length: 
    my_protein = ampal.load_pdb(f'{prot_name}.pdb')
    chain_A_length = len(my_protein[0].sequence)
    split_dict = {}
    for mut_name, full_seq in mutated_dict.items():
        chain_A = full_seq[:chain_A_length]
        chain_B = full_seq[chain_A_length:]

        split_dict[mut_name] = {
            "chain_A": chain_A,
            "chain_B": chain_B
        }


    # Save the mutated sequences to a JSON file
    with open(f'{prot_name}_doublechain_mut_all.json', 'w') as f:
        json.dump(split_dict, f, indent=4)
    print(f"Mutated combined chain sequences saved to '{prot_name}_doublechain_mut_all.json'.")
    
# example usage

if __name__ == "__main__":
    prot_name = "4dri_pose5_corrected" # adjust the name of protein
    pocket_size = 10.0
    cutoff = 0.135
    do_seq_design_pipeline(prot_name, pocket_size, cutoff)