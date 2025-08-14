#Function which runs boltz, must specify YAML folder. 
#Based on a script by Handing Wang, Wells-Wood-Research. 
#Updated by student B144319
import os 


def run_boltz_batch_prediction_sequentially(input_yaml_folder: str):
    """
    Runs Boltz prediction sequentially for all YAML files within a specified folder.

    Args:
        input_yaml_folder (str): The path to the folder containing batch YAML files.

    """
    current_dir = os.getcwd()
    # Get a list of all YAML files in the input folder
    yaml_files = [f for f in os.listdir(input_yaml_folder) if f.endswith('.yaml')]

    for yaml_file in yaml_files:
        yaml_path = os.path.join(input_yaml_folder, yaml_file)
        # check if the boltz results already exist, if already finished, skip
        output_pdb_path = os.path.join(current_dir, "boltz_results",yaml_file.replace('.yaml', '.pdb'))
        if os.path.exists(output_pdb_path):
            print(f"Skipping {yaml_file}, already processed.")
            continue

        # boltz command in shell
        boltz_command = f"boltz predict {yaml_path}  --use_msa_server --output_format pdb"
        
        os.system(boltz_command)
        
 
run_boltz_batch_prediction_sequentially("")

