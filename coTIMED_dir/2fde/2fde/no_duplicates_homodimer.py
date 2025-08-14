import os 
import json 

input_dir = "/home/s_alya/research_project/chosen_prots_timed/2fde/2fde"
prot_name = "2fde"
unique_seqs={}
seen = set()

for file in os.listdir(input_dir):
    if file.endswith('_mut_all.json'):
        file_path = os.path.join(input_dir, file)
        
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        for file_name, sequence in data.items():
            if sequence not in seen:
                unique_seqs[file_name] = sequence
                seen.add(sequence)
            else:
                print(f'Duplicate removed in {file}: {file_name}')


output_path = os.path.join(input_dir, f'unique_sequences_{prot_name}.json')
with open(output_path, 'w') as f:
    json.dump(unique_seqs, f, indent=4)

print(f'Finished processing files for {prot_name}.')
