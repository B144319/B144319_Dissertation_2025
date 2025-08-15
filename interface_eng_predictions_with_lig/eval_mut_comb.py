__version__ = "0.1.5"
# updating 2025-07-07

import json
import pandas as pd 
import itertools
from collections import defaultdict
import random
import os
import numpy as np

# Import PyRosetta modules
import pyrosetta
from pyrosetta import *
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, InitializeFromCommandline
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueIndexSelector,
    NeighborhoodResidueSelector
)
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent,
    OperateOnResidueSubset,
    PreventRepackingRLT,
    RestrictToRepackingRLT
)

'''

updating 2025-05-29
The script reads a csv file containing the probablities of each aminio acid at each position,
filtering out positions with relatively high probabilities for any amino acid.
Author: HD Wang
Arguments:
    - `redesign_path`: Path to the CSV file containing amino acid probabilities. From cotiemed result.
    - `pdb_path`: Path to the PDB file of the wild-type structure.
    - `pose`: PyRosetta pose object of the homodimer.
Returns:
    - `filtered_aas`: Dictionary of filtered amino acids with positions as keys and lists of substitutable amino acids as values.
    - `obvious_mutations`: Dictionary of mutations with ddG > 5 kcal/mol.

'''


def setup_relax_mover(scorefxn, mut_position = None):
    """
    Set up the movers for relaxation: packer and minimizer
    
    Args:
        scorefxn: PyRosetta score function
        mut_position: Residue position to mutate (pose numbering)
    
    Returns:
        tuple: (packer mover, minimizer mover)
    """
    

    # Create deterministic task factory
    tf = TaskFactory()

    # These are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Restrict packing to only mutation sites and neighbors
    if mut_position:
        mut_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(",".join(map(str, mut_position)))
    
        neighbor_selector = NeighborhoodResidueSelector(mut_selector, distance=6.0, include_focus_in_subset=True)
        
        # Allow repacking only in selected residues
        tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), neighbor_selector))
        tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), ~neighbor_selector))
    
    else: # For WT pre-relax: allow full repacking
        # If no mutation, just repack all
        tf.push_back(RestrictToRepacking())


    # Configure deterministic packer
    packer = PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.nloop(6)  # Increased packing cycles for convergence

    # Configure tight minimization
    move_map = MoveMap()
    move_map.set_bb(True)  # Allow backbone movement
    move_map.set_chi(True)
    move_map.set_jump(False)  # No jumps in homodimer
    
    min_mover = MinMover()
    min_mover.score_function(scorefxn)
    min_mover.movemap(move_map)
    min_mover.min_type("lbfgs_armijo_nonmonotone")  # More stable than dfpmin
    min_mover.tolerance(0.001)  # Tighter convergence
    min_mover.max_iter(1000)    # Increased iterations
    
    return packer, min_mover


# udpating the cal of wt_ddG_interface
def pre_relax_wt(pose):


    """Pre-relax the wild-type structure once
    Args:
        pose: PyRosetta pose object of the wild-type structure
    Returns:
        pose: Relaxed PyRosetta pose object of the wild-type structure
    
    2025-05-29, Updated
    """
    scorefxn = pyrosetta.get_score_function()
    wt_pose = pose.clone()
    
    # Full repack and minimization for initial relaxation
    packer, min_mover = setup_relax_mover(scorefxn)
    
    # Add constraints to maintain overall structure
    add_constraints = pyrosetta.rosetta.protocols.constraint_generator.AddConstraints()
    
    add_constraints.apply(wt_pose)
    
    # Apply 3 cycles of pack-minimize
    for _ in range(3):
        packer.apply(wt_pose)
        min_mover.apply(wt_pose)
    
    return wt_pose


def cal_ddg_interface(pre_releaxed_wt, mutation, position, relax=True):
    '''
    Calculate ddG for a symmetric homodimer with mutations at the same position in both chains.
    
    Args:
        pose: PyRosetta pose object of the homodimer
        mutation: Single-letter amino acid code to mutate to
        position: PDB numbering position to mutate
        relax: Whether to perform side-chain repacking and minimization
    
    Returns:
        tuple: (ddG, wild-type dG, mutant dG)
    
    2025-04-30, Updated
    '''   
    # Store original pose
    wt_pose = pre_releaxed_wt.clone()
    

    try:
                

        # Convert PDB to pose numbering
        resA = pre_releaxed_wt.pdb_info().pdb2pose('A', position)
        resB = pre_releaxed_wt.pdb_info().pdb2pose('B', position)
        if resA == 0 or resB == 0:
            raise ValueError(f"Residue {position} not found in both chains!")

        # Create mutant pose
        mutant_pose = pre_releaxed_wt.clone()
        mutate_residue(mutant_pose, resA, mutation)
        mutate_residue(mutant_pose, resB, mutation)

        # Use Cartesian score function 
        #scorefxn = pyrosetta.create_score_function("ref2015_cart")
        #scorefxn.set_weight(pyrosetta.rosetta.core.scoring.cart_bonded, 1.0)  # Add cart_bonded term


        # Use standard score function for stability
        scorefxn = pyrosetta.get_score_function()  # Instead of ref2015_cart， here ues ref2015
        
        if relax:
            # Setup movers with focused repacking
            packer, min_mover = setup_relax_mover(scorefxn, [resA, resB])

            # Add constraints to maintain structure
            add_constraints = pyrosetta.rosetta.protocols.constraint_generator.AddConstraints()
            
            add_constraints.apply(mutant_pose)

            # Apply identical relaxation to both states
           
            packer.apply(mutant_pose)
            for _ in range(3):  # Multiple minimization cycles
                min_mover.apply(mutant_pose)

        # Energy calculation with consistent settings
        ia = InterfaceAnalyzerMover("A_B")
        ia.set_scorefunction(scorefxn)
        ia.set_pack_separated(True)
        

        # Calculate wild-type energ
        wt_ia = ia.clone()
        wt_ia.apply(wt_pose)
        wt_dG = wt_ia.get_interface_dG()

        # Calculate mutant energy
        mut_ia = ia.clone()
        mut_ia.apply(mutant_pose)
        mutant_dG = mut_ia.get_interface_dG()

        ddG = mutant_dG - wt_dG
        print(f"ΔΔG: {ddG:.2f} REU (WT: {wt_dG:.2f}, Mut: {mutant_dG:.2f})")


        # get the mutant sequence
        mutant_sequence = mutant_pose.sequence()
        
        return ddG, wt_dG, mutant_dG,mutant_sequence
    except Exception as e:
        raise RuntimeError(f"ddG calculation failed: {str(e)}") 
    

def get_substitution_aa(df, pose):
    """Extract substitution amino acids from a row of the DataFrame.


    Args:
        df (pd.DataFrame): DataFrame containing amino acid probabilities. from cotiemed result
        the format is strictly as it only process 1 protien each time

        updated 20250609: don't forget the chainB residue would add the length of chainA

        pose (pyrosetta.Pose): PyRosetta pose object to get wild-type amino acids.
    Returns:
        dict: A dictionary with position as keys and lists of substitutable amino acids as values.
    2025-05-29, Updated
    
    """

    filtered_aas = {}
    for idx, row in df.iterrows():


        amino_idx = int(row['idx'])
        print(f"Processing position {amino_idx}...")
        
    
        wildtype_aa = pose.residue(amino_idx).name3() # Get the wildtype amino acid at this position

        wt_id = wildtype_aa + str(amino_idx) # Create a unique identifier 
        # get subsitution aa with probability > 0.1

        filtered_aa = {}
        substitution_aa = []
        for col in df.columns[2:]:  # Skip the first 2 columns which is the name position

            if row[col] >= 0.1:  # choose high probability amino acids
                substitution_aa.append(col)


        filtered_aa[wt_id] = substitution_aa
        print(f"Position {amino_idx} can be substituted with: {filtered_aa}")

        filtered_aas.update(filtered_aa)

    return filtered_aas


def do_filtered_aas(filtered_aas, pre_relaxed_wt):
    """Filter amino acids based on specific criteria.
    
    Args:
        filtered_aas (dict): Dictionary of filtered amino acids.
        e.g. {'MET1': ['S', 'T'], 'LYS2': ['Q', 'R']}
        
    Returns:
        dict: Filtered amino acids dictionary.
    """
    obvious_mutations = {}  # initialize the dictionary
    auto_mutations = {}  # initialize the dictionary for auto mutations
    medium_mutations = {}  # initialize the dictionary for medium mutations
    for position, aas in filtered_aas.items():
        print(f"\nProcessing position {position} with substitutions: {aas}")
        
        if not aas:
            continue  # skip if no substitutions available
            
        try:
            # Extract residue number (works for both 'MET1' and 'A1' formats)
            res_num = int(''.join(filter(str.isdigit, position)))
            
            # Initialize position entry if not exists
            if position not in obvious_mutations:
                obvious_mutations[position] = {}
            if position not in auto_mutations:
                auto_mutations[position] = {}   
            if position not in medium_mutations:
                medium_mutations[position] = {}

            for mutant_aa in aas:
                try:
                    ddG, wt_dG, mutant_dG,mutant_sequence = cal_ddg_interface(
                        pre_relaxed_wt, 
                        mutant_aa, 
                        res_num,  # pass the integer residue number
                        relax=True
                    )
                    
                    print(f"{position}{mutant_aa}: ddG = {ddG:.2f} kJ/mol")

                    # store possible mutations into auto_mutations no matter what the ddG is
                    auto_mutations[position][mutant_aa] = {
                        "ddG": float(f"{ddG:.2f}"),
                        "wt_dG": float(f"{wt_dG:.2f}"),
                        "mutant_dG": float(f"{mutant_dG:.2f}"),
                        "mutant_sequence": mutant_sequence
                    }
                    if ddG > 2:  # Filter for ddG > 2
                        medium_mutations[position] = {
                            "mutant_aa": mutant_aa,
                            "ddG": ddG,
                            "wt_dG": wt_dG,
                            "mutant_dG": mutant_dG,
                            "mutant_sequence": mutant_sequence
                        }



                    # Only store obvious mutations with ddG > 5

                    if ddG > 5:  # Changed from 5 to 0 based on your code
                        obvious_mutations[position][mutant_aa] = {
                            "ddG": float(f"{ddG:.2f}"),
                            "wt_dG": float(f"{wt_dG:.2f}"),
                            "mutant_dG": float(f"{mutant_dG:.2f}"),
                            "mutant_sequence": mutant_sequence
                        }
                        
                    

                except Exception as e:
                    print(f"Error processing {position}{mutant_aa}: {str(e)}")
                    continue
                    
        except Exception as e:
            print(f"Error processing position {position}: {str(e)}")
            continue
            
    return auto_mutations, medium_mutations, obvious_mutations



def do_mut_pipeline(pdb_path: str, redesign_path: str, output_prefix: str = "demo",cutoff: float = 2.0):
    """
    Main function to run the mutation pipeline.
    
    Args:
        pdb_path (str): Path to the PDB file of the wild-type structure.
        redesign_path (str): Path to the CSV file containing amino acid probabilities.
        output_prefix (str): Prefix for output files.
        
    Returns:
        json file: includes the predicted mutations with ddG > 2.
                   abovious mutations with ddG > 5.
                   and all the mutations with ddG > 0, named auto_mutations.json.

    2025-05-29, Updated
    """
    # Initialize PyRosetta
    pyrosetta.init()

    # Create output directory if it doesn't exist
    os.makedirs(output_prefix, exist_ok=True)

    pose = pyrosetta.pose_from_pdb(pdb_path)
    pre_relaxed_wt = pre_relax_wt(pose)
    wt_sequence = pre_relaxed_wt.sequence()
    print(f"Wild-type sequence: {wt_sequence}")

    df = pd.read_csv(redesign_path)
    filtered_aas = get_substitution_aa(df, pre_relaxed_wt)
    print(filtered_aas)
    try:
        auto_mutations, medium_mutations, obvious_mutations = do_filtered_aas(filtered_aas, pre_relaxed_wt)

        # Write the obvious mutations to a JSON file        


        with open(f"{output_prefix}/{output_prefix}_obvious_mutations.json", "w") as f:
            json.dump(obvious_mutations, f, indent=4)
        # Write the auto mutations to a JSON file
        with open(f"{output_prefix}/{output_prefix}_auto_mutations.json", "w") as f:
            json.dump(auto_mutations, f, indent=4)

        # filter the mutations with ddG >2, and throw the empty mutations:
        predicted_mutations = {}
        for position, mutations in obvious_mutations.items():
            for mutant_aa, data in mutations.items():

                # if there has 2 or more mutations at the same position, we need to include all of them
                if position not in predicted_mutations:
                    predicted_mutations[position] = {}
                # if the ddG is greater than cutoff, we will keep otherwise we will not keep it
                if data["ddG"] > cutoff:
                    # if the mutant_aa is not in the predicted_mutations, we will add it;if the mutant_aa is already in the predicted_mutations, we will not overwrite it
                    if mutant_aa not in predicted_mutations[position]:
                        predicted_mutations[position][mutant_aa] = {
                            
                            "ddG": data["ddG"],
                            "wt_dG": data["wt_dG"],
                            "mutant_dG": data["mutant_dG"],
                            "mutant_sequence": data["mutant_sequence"]
                        }
        # remove empty mutations
        predicted_mutations = {k: v for k, v in predicted_mutations.items() if v}

        with open(f"{output_prefix}/{output_prefix}_predicted_mutations.json", "w") as f:
            json.dump(predicted_mutations, f, indent=4)
        
        return predicted_mutations
    except Exception as e:
        print(f"Error in mutation pipeline: {str(e)}")

        # return an empty dict if any error occurs
        return {}

###########################
######Section2 : get score from the mutation interface

def get_interface_scores(pdb_file):

    
    """
    Calculate interface scores for a given pose.
    
    Args:
        pose (pyrosetta.Pose): PyRosetta pose object.
        chain1 (str): First chain identifier (default 'A').
        chain2 (str): Second chain identifier (default 'B').
        
    Returns:
        dict: Dictionary containing ptm, iptm, and per_chain_ptm scores.
    """
        # load pose
    pose = pyrosetta.pose_from_pdb(pdb_file)
    
    # call the relax function to pre-relax the pose
    pre_relaxed_wt = pre_relax_wt(pose)

    # get the interface ddG



    # analyze interface statistics
    ia = InterfaceAnalyzerMover()
    ia.set_interface("A_B")
    scorefxn = pyrosetta.get_fa_scorefxn()
    ia.set_scorefunction(scorefxn)
    ia.set_compute_packstat(True)
    ia.set_compute_interface_energy(True)
    ia.set_calc_dSASA(True)
    ia.set_calc_hbond_sasaE(True)
    ia.set_compute_interface_sc(True)
    ia.set_pack_separated(True)
    ia.apply(pre_relaxed_wt)


  

    # retrieve statistics
    interfacescore = ia.get_all_data()
    interface_sc = interfacescore.sc_value # shape complementarity
    interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
    interface_dG = ia.get_interface_dG() # interface dG
    interface_dSASA = ia.get_interface_delta_sasa() # interface dSASA (interface surface area)
    interface_packstat = ia.get_interface_packstat() # interface pack stat score
    #  ratio of dG/dSASA*100 is most useful for interface energy normalisation
    interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100 # ratio of dG/dSASA (normalised energy for interface area size)

    return {
        'dG': interface_dG,
        'sc': interface_sc,
        'interface_hbonds': interface_interface_hbonds,
        'interface_packstat': interface_packstat,
        'dSASA': interface_dSASA,
        'dG_SASA_ratio': interface_dG_SASA_ratio,

    }



###############
##### section3: filter the mutations and generate the comb mutations
def get_single_chain_sequence(pdb_path: str) -> str:
    """
    Get the sequence of the first chain from a PDB file.
    
    Args:
        pdb_path (str): Path to the PDB file.
        
    Returns:
        str: Sequence of the first chain.
    """
    # Initialize PyRosetta
    pyrosetta.init()
    pose = pyrosetta.pose_from_pdb(pdb_path)
    return pose.chain_sequence(1)  # Get the sequence of the first chain


def filter_mutations(ddg_mutation: dict, min_ddg: float = 2, max_ddg: float = 10) -> list:
    mutation_options = []

    for position, sub_dict in ddg_mutation.items():
        resi = int(position[3:])
        for mutant_aa in sub_dict.keys():
            try:
                mut_data = sub_dict[mutant_aa]
                ddG = mut_data.get('ddG')
                if min_ddg < ddG < max_ddg:
                    mutation_options.append((resi, mutant_aa))  # <-- tuple, not list of tuple
            except (ValueError, AttributeError) as e:
                print(f"Skipping {position}: {str(e)}")
                continue
    
    return mutation_options

def generate_combinations(mutation_options: list, min_mutations: int, max_mutations: int, num_combinations: int = 1000):
    """
    Generate random combinations of mutations without repeating positions in each combination.

    Args:
        mutation_options (list): List of lists of mutations (each mutation is a tuple like (position, amino_acid))
        num_combinations (int): Number of random combinations to generate
        min_mutations (int): Minimum number of mutations per combination
        max_mutations (int): Maximum number of mutations per combination

    Returns:
        list: List of randomly generated mutation combinations
        format example:
        [
            [(4, 'V'), (5, 'I')],
            [(4, 'V'), (6, 'A')],
            ...
        ]
    """
    # Group mutations by position
    position_dict = defaultdict(list)
    for sublist in mutation_options:
        for mutation in sublist:
            position_dict[mutation[0]].append(mutation)

    # Get unique positions
    positions = list(position_dict.keys())
    
    # Ensure we don't ask for more mutations than unique positions
    max_mutations = min(max_mutations, len(positions))
    min_mutations = min(min_mutations, max_mutations)
    
    all_valid_combinations = []

    for r in range(min_mutations, max_mutations + 1):
        for pos_combo in itertools.combinations(positions, r):
            # For each position in pos_combo, pick one of its mutation options
            mutation_sets = [position_dict[pos] for pos in pos_combo]
            for mutation_choice in itertools.product(*mutation_sets):
                all_valid_combinations.append(mutation_choice)

    # Sample if needed as we probably don't need all combinations we can test them all 

    if len(all_valid_combinations) > num_combinations:
        return random.sample(all_valid_combinations, num_combinations)
    else:
        return all_valid_combinations


#step3, for all the combinations, generate the fasta file - edited for heterodimers (non-identical chain A and B mutations)

def generate_fasta(combinations: list, wt_sequences: dict, save_path: str, write_comb_info: bool = True, write_fasta: bool = True):
    """
    Generate FASTA files and mutation combination metadata for heterodimers.

    Args:
        combinations (list): List of mutation combinations, each as [(chain, pos, aa), ...]
        wt_sequences (dict): Dict of wild-type sequences, e.g., {'A': '...', 'B': '...'}
        save_path (str): Path to save FASTA and combination info
        write_comb_info (bool): Whether to save mutation metadata (default: True)
        write_fasta (bool): Whether to write FASTA files (default: True)
    """
    import os
    import json

    if (write_fasta or write_comb_info) and not save_path:
        raise ValueError("save_path must be provided when write_fasta or write_comb_info is True")

    os.makedirs(save_path, exist_ok=True)

    comb_info = {}

    for i, combo in enumerate(combinations):
        # Start from wild-type sequences
        mutated_sequences = {chain: list(seq) for chain, seq in wt_sequences.items()}

        # Apply mutations to appropriate chains
        for chain, pos, aa in combo:
            if chain in mutated_sequences:
                index = pos - 1  # 1-based to 0-based
                mutated_sequences[chain][index] = aa
            else:
                raise ValueError(f"Unknown chain ID '{chain}' in mutation combination.")

        # Join sequences back to strings
        mutated_sequences = {chain: ''.join(seq) for chain, seq in mutated_sequences.items()}

        comb_info[f"comb_{i + 1}"] = {
            "mutations": combo,
            "heterodimer_sequence": mutated_sequences
        }

        if write_fasta:
            fasta_path = os.path.join(save_path, f"comb_{i + 1}.fasta")
            with open(fasta_path, "w") as f:
                for chain_id in sorted(mutated_sequences):
                    f.write(f">protein|name=chain_{chain_id}\n{mutated_sequences[chain_id]}\n")

    if write_comb_info:
        with open(os.path.join(save_path, "comb_info.json"), "w") as json_file:
            json.dump(comb_info, json_file, indent=4)

    print(f"Successfully processed {len(combinations)} combinations.")



def parse_prot_mut_plddt(wt_pdb, design_pdb, plddt_file):
    """
    Parse the plddt values for the mutated residues in the design pdb.
    design_pdb better remove the ligand, as ligand atom would missed in boltz result causing 
    pyrosequence error.

    Args:
        wt_pdb (str): Path to the wild-type PDB file.
        design_pdb (str): Path to the redesign PDB file.
        plddt_file (str): Path to the plddt.npz file.
    """
    data = np.load(plddt_file)
    plddt = data['plddt.npy']
    
    wt_seq = get_single_chain_sequence(wt_pdb)
    design_seq = get_single_chain_sequence(design_pdb)

    mutat_idx = []
    for i, (wt_res, design_res) in enumerate(zip(wt_seq, design_seq)):
        if wt_res != design_res:
            mutat_idx.append(i)

    plddt_mutat_values = []
    plddt_mutat_info = {}
    for i in mutat_idx:
        plddt_mutat_values.append(plddt[i])
        plddt_mutat_info[i] = plddt[i]

    mean_mut_plddt = np.mean(plddt_mutat_values)
    
    return mean_mut_plddt, plddt_mutat_info


# Example usage

if __name__ == "__main__":
    # Example paths, replace with your actual paths
    pdb_path = '/home/hdwang/protein_sensor_hd/1flm/pdb/1flm_rmlig.pdb'

    redesign_path= "/home/hdwang/protein_sensor_hd/1flm/timed_design/separate_interface/1flm_HCY4_predict_prob.csv"

    output_prefix = "1flm_hcy4"

    # Run the mutation pipeline
    do_mut_pipeline(pdb_path, redesign_path, output_prefix)
    
    print("Mutation pipeline completed successfully.")