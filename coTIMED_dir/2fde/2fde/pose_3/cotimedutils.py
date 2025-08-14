__version__ = "0.1.2"

import pandas as pd
import json
import pathlib
import os
from itertools import product
import ampal
from pathlib import Path
import os
from collections import OrderedDict
from collections import defaultdict
from Bio import SeqIO
import ampal.geometry
import numpy as np
import json

def residue_number_indexing(assembly: ampal.Assembly):
    count = 0
    residue_number = []

    for chain in assembly:
        for residue in chain:
            count += 1
            residue.id = count
            residue_number.append(int(residue.id))

def get_pocket_residues(cutoff) -> dict:
    """
    Extracts pocket residues from a PDB file and returns them as a dictionary.
    
    Args:
        all pdb files under current directory
        cutoff (float): The distance cutoff to consider residues as surrounding the ligand.

    Returns:
        dict: A dictionary where keys are chain identifiers and values are lists of residue indices
        surrounding the ligand.
    """

    current_directory = os.getcwd()

    files = [f for f in os.listdir(current_directory) if f.endswith(".pdb")]
    ligand_collection = []

    with open("/home/s_alya/research_project/ligand.csv", "r") as ligand_file:
        for row in ligand_file:
            ligand_collection.append(row)

    ligand_collection = [x.strip("\n") for x in ligand_collection]


    ligand_surrounding_res_dict = defaultdict(list)

    all_compressed_chains = []
    all_compressed_residues = []
    countlig = []
    countfile = 0
    backbone_atoms = ["C", "N", "O", "CA"]

    for file in files:
        protein = ampal.load_pdb(file)
        _ = residue_number_indexing(protein)
        ligands = protein.get_ligands()
        
        for ligand in ligands:
            print(ligand)
            if (
                ligand.mol_code in ligand_collection
            ):  # Assuming ligand_collection is defined somewhere
                ligand_atoms = ligand.atoms.keys()
                residues_close_to_ligand = []
                chains_close_to_ligand = []
                # Use defaultdict to store residues grouped by chains
                chain_residue_map = defaultdict(set)

                # Iterate over ligand atoms and find close atoms in protein
                for atom_name in ligand_atoms:
                    for atom in protein.get_atoms(ligands=False):
                        if atom.res_label not in backbone_atoms:
                            if ampal.geometry.distance(ligand[atom_name], atom) <= cutoff:
                                residues_close_to_ligand.append(atom.parent.id)
                                chains_close_to_ligand.append(atom.parent.parent.id)

                # Group residues by chains using defaultdict

                for residue, chain in zip(residues_close_to_ligand, chains_close_to_ligand):
                    base_filename = os.path.splitext(file)[0]
                    name = f'{base_filename}_{chain}'
                    chain_residue_map[name].add(int(residue))

                # Convert defaultdict values to lists
                compressed_chains = list(chain_residue_map.keys())
                compressed_residues = [
                    list(residues) for residues in chain_residue_map.values()
                ]

                # Append to the overall results
                all_compressed_chains.extend(compressed_chains)
                all_compressed_residues.extend(compressed_residues)


    for residue, chain in zip(all_compressed_residues, all_compressed_chains):
        ligand_surrounding_res_dict[chain] = residue

    filename = 'pocket_residues.json'
    with open(filename, 'w') as f:
        json.dump(ligand_surrounding_res_dict, f, indent=4)

    return ligand_surrounding_res_dict



def process_data(probabilities_file: str, label_map_file: str, output_file: str) -> pd.DataFrame:
    """
    Processes and combines probabilities data with label map data.

    Args:
        probabilities_file (str): Path to the CSV file containing probabilities.
        label_map_file (str): Path to the text file containing label mappings.
        output_file (str): Desired name for the output CSV file.

    Returns:
        pd.DataFrame: The combined and processed DataFrame.
    """
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'

    # Read probabilities data
    df_1 = pd.read_csv(probabilities_file, header=None, names=[i for i in alphabet])

    # Read label map data
    df_2 = pd.read_csv(label_map_file, header=None, names=['pdbinfo', 'chain', 'idx', 'AAs'])
    df_2['pdbinfo_chain'] = df_2['pdbinfo'] + '_' + df_2['chain']

    # Combine the two dataframes
    df_3 = pd.concat([df_2['pdbinfo_chain'], df_2['idx'], df_1], axis=1)

    # Save the combined dataframe to a CSV file
    df_3.to_csv(output_file, index=False)

    print(f"Data successfully processed and saved to '{output_file}'")
    return df_3


def get_amino_acid(current_pdb_chain: str, residue_idx: int, df: pd.DataFrame, cutoff) -> list:
    """
    Identifies amino acids with probabilities above a given cutoff for a specific
    protein-chain and residue index.

    Args:
        current_pdb_chain (str): A string in the format 'pdbid_chain' (e.g., '3cyl_A').
        residue_idx (int): The residue index.
        df (pd.DataFrame): The DataFrame containing amino acid probabilities and residue information.
        cutoff (float, optional): The probability cutoff. Defaults to 0.1.

    Returns:
        list: A list of amino acids meeting the cutoff criteria.
    """
    
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'

    # Efficiently filter the DataFrame
    selected_row = df[(df['pdbinfo_chain'] == current_pdb_chain) & (df['idx'] == residue_idx)]

    selected_residues = []
    if not selected_row.empty:
        # Extract probabilities for the selected row
        probabilities = selected_row[list(alphabet)].iloc[0]
        for aa in alphabet:
            if probabilities[aa] >= cutoff:
                selected_residues.append(aa)

    return selected_residues


def generate_fasta(protein_path):
    '''
    For a given protein pdb file, this function generates a dictionary
    e.g. {1: 'A', 2: 'C', 3: 'D', ...}
    
    '''


    # Load the PDB file
    my_structure = ampal.load_pdb(protein_path)

    # get fasta sequence
        
    sequences = my_structure.sequences
    
    # extract the sequence from the list and convert it to a string
    fasta_seq = ''.join(sequences)
    
    # create a dictionary to store the residue number and the corresponding three letter code
    residue_dict = {}
    for i in range(0,len(fasta_seq)):
        residue_dict[i+1] = fasta_seq[i]
    
    return residue_dict


def generate_combination(prot_name: str, predicted_json_path: str) -> list:

    ''' 
        In cotimed, which mutate all residues the idx is from 1 to total number of residues,
        so we need to relabel the residue number, here I generate new dict to store all residues

        Arguments:
        prot_name: the name of the protein, e.g. '1flm_hcy12'
        predicted_json_path: the path to the predicted JSON file containing mutation predictions
        e.g. 

        residue_dict: a dictionary containing residue numbers and their corresponding amino acids
        e.g. {1: 'A', 2: 'C', 3: 'D', ...}
        This function generates all possible combinations of mutations for the specified protein.

        returns:
        permutation_matrix: a list of dictionaries, each dictionary represents a combination of mutations
        e.g. [{'1': 'A', '2': 'C', '3': 'D'}, {'1': 'A', '2': 'C', '3': 'E'}, ...]
    '''
    
    print(f'Processing {prot_name}')
    mutate_dict = {}
    
    with open(predicted_json_path, 'r') as f:
        all_predict = json.load(f)


    for key in all_predict:

        if prot_name in key:

            ''' no matter how many chains, we introduce all residues'''

            mutate_dict.update({k: v for k, v in all_predict[key].items() if v != []})
  
    # create permutation matrix to store the residue number that need to be mutated
    permutation_matrix = []

    # 1. sort the keys in the mutate_dict, and get the order of the residue position
    sorted_keys = sorted(mutate_dict.keys()) # create a list of position keys

    candidates_list = [mutate_dict[key] for key in sorted_keys]
 
    # 2.generate all possible combinations
    all_combinations = list(product(*candidates_list))


    # 3. add the residue number to the permutation matrix

    for i in all_combinations:
        temp = {}
        for j in range(len(i)):

            temp[sorted_keys[j]] = i[j]

        permutation_matrix.append(temp)

    return permutation_matrix

# for heterodimer combine both chains into a single chain
def get_single_chain(pdb_path: str, ) -> str:
    """
    Extracts the amino acid sequence from a PDB file and returns it as a string.
    Args:
        pdb_path (str): Path to the PDB file.
    Returns:
        str: Amino acid sequence extracted from the PDB file.
    """
    # Load the PDB file
    my_structure = ampal.load_pdb(pdb_path)

    # get fasta sequence
    sequences = my_structure.sequences
    
    # extract the sequence from the list and convert it to a string
    fasta_seq = sequences[0]
    
    return fasta_seq


# mutate all residues, iterate the dict

def all_on_1chain(all_predict: dict, sigle_chain: str,prot_new_name: str)-> dict:
    
    '''
    As a homodimer, we need to mutate all residues in the single chain, when we synthesize the protein.

    Arguments:
            all_predict: dict, a json file that contains all residues that need to be mutated
            sigle_chain: str, the fasta sequence of the sigle chain generted from the pdb file use ampal

            prot_new_name: str, the new name of the protein, e.g. '1flm_hcy12'
    
    return: mutated_dict_list: list, a list of dict, each dict contains the new fasta sequence, save to a json file

    revised by H.D.Wang 2025.04.12
    
    '''
   
    num = 0
    mutated_dict_all = {}
    for i in all_predict:
        
        num += 1
        redesign_prot_name  = prot_new_name + '_' + str(num)
    
        # use dict to store new fasta sequence
        
        fasta_seq_copy = sigle_chain # start with the original single chain sequence
        for key in i:
            mutated_dict_iter = {}
            # get the residue number
            residue_num = int(key)
            # get the new residue
            residue_new = i[key]


            # mutate the residue in the sigle chain
            if residue_num <= len(sigle_chain):
                fasta_seq_copy = fasta_seq_copy[:residue_num-1] + residue_new + fasta_seq_copy[residue_num:]
            
            else:
                new_idx = residue_num - len(sigle_chain)
                fasta_seq_copy = fasta_seq_copy[:new_idx-1] + residue_new + fasta_seq_copy[new_idx:]
        
        # store all new fasta sequence in the dict
        mutated_dict_iter[redesign_prot_name] = fasta_seq_copy
        # store the dict in the list
        mutated_dict_all.update(mutated_dict_iter)

    return mutated_dict_all

# section2: writing the sequence from onchain_mut_all.json to fasta files.

def write_fasta(seq_list,filename,save_path):

    path_2_save = f'{save_path}/{filename}'
    with open(path_2_save, 'w') as f:
        header_prefix = '>protein|name=chain_'
        header_list = [header_prefix + chr(65 + i) for i in range(len(seq_list))]

        for seq,header in zip(seq_list,header_list):
            # add header
            f.write(header + '\n')
            # write sequence
            f.write(seq + '\n')


def write_fasta_from_dict(fasta_matrix, save_folder: str):
    '''
    
    input: fasta_matrix_folder: str, the path of the fasta sequence matrix
           save_folder: str, the path of the folder to save the fasta sequence
    output: save the fasta sequence to the save_folder
    description:
    1. read the fasta sequence matrix from result generated by all_on_1chain
    2. extract the key as filename
    3. extract the value as sequence
    4. save the sequence to a fasta file, as 2fnu is a dimer, so save the sequence twice
    5. check if the save folder exists if not create it
    6. write the fasta file
    
    revised by H.D.Wang 2025.04.12
    '''

        # open the fasta sequence matrix

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    
    for filename, sequence in fasta_matrix.items():
        # Save as <filename>.fasta
         
        write_fasta([sequence,sequence], filename,save_folder)  # Wrap sequence in a list for write_fasta
        
        




def do_seq_design_pipeline(protein_name: str,pdb_path: str,cutoff) -> list:

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

    pdb_name = protein_name.split('.')[0] # e.g. '3cyl_TES_A_124' in the form of pdbid_lig_ligchain_ligid
    probabilities_path = '/home/s_alya/research_project/chosen_prots_timed/2fde/2fde/pose_3/CoTIMED.csv'
    label_map_path = 'chosen_prots_timed/2fde/2fde/pose_3/datasetmap.txt'
    output_csv_name = f'{pdb_name}_predict_prob.csv'

    # Process data and get the DataFrame
    df_prob = process_data(probabilities_path, label_map_path, output_csv_name)

    with open(f'{pdb_name}_pocket_residues.json', 'w') as f:
        json.dump(get_pocket_residues(cutoff), f, indent=4)
        print(f"Processed probabilities data saved to '{output_csv_name}'.")

    # --- Section 2: Process the redesign sequence and JSON file from the selected sequence ---


    # Load redesign sequence information from a JSON file

    redesign_sequence = get_pocket_residues(cutoff)
    current_directory = os.getcwd()
    
    predict_seq = {}
    for key_full_name in redesign_sequence.keys():

        predict_seq[key_full_name] = {}  # Initialize a sub-dictionary for each key_full_name
        for item_idx in redesign_sequence[key_full_name]:
            # Ensure item_idx is an integer for lookup
            predict_seq[key_full_name][item_idx] = get_amino_acid(key_full_name, int(item_idx), df_prob, cutoff)

    # Save the predicted sequence to a JSON file
    output_json_name = f'{pdb_name}_redesign_pocket_sequence.json'
    json_output_path = os.path.join(current_directory, output_json_name)

    with open(json_output_path, 'w') as f:
        json.dump(predict_seq, f, indent=4)
    print(f"Predicted sequence data saved to '{output_json_name}'.")

    # --- Section 3: Generate the mutatation combinations based on the predicted JSON file ---
    # Generate the mutation combinations
    redesign_sequence = generate_fasta(pdb_path)
    mutation_combinations = generate_combination(pdb_name, json_output_path)
    
    # as a homodimer, we need to mutate all residues in the single chain, when we synthesize the protein.
    single_chain = get_single_chain(pdb_path)
    mutated_dict_list = all_on_1chain(mutation_combinations, single_chain, pdb_name)

    # Save the mutated sequences to a JSON file
    with open(f'{pdb_name}_onechain_mut_all.json', 'w') as f:
        json.dump(mutated_dict_list, f, indent=4)
    print(f"Mutated sequences saved to '{pdb_name}_onechain_mut_all.json'.")
    
# example usage
if __name__ == "__main__":

    prot_name = "/home/s_alya/research_project/chosen_prots_timed/2fde/2fde/pose_3/2fde_pose3_corrected" # adjust the name of protein
    pocket_size = 10.0 #increased pocket size to 10. 
    cutoff = 0.115 #changed cutoff from 0.1 
    do_seq_design_pipeline(prot_name, pocket_size, cutoff)

