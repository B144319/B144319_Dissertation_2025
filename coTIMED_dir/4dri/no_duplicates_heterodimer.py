import os
import json

input_dir = "/home/s_alya/research_project/chosen_prots_timed/4dri"
prot_name = "4dri"

unique_seqs = {}
seen_pairs = set()

for file in os.listdir(input_dir):
    if file.endswith('_mut_all.json'):
        file_path = os.path.join(input_dir, file)

        with open(file_path, 'r') as f:
            data = json.load(f)

        for mut_name, chains in data.items():
            chain_a = chains.get("chain_A", "")
            chain_b = chains.get("chain_B", "")

            # Treat (chain_A, chain_B) as an ordered tuple
            pair_key = (chain_a, chain_b)

            if pair_key not in seen_pairs:
                unique_seqs[mut_name] = chains
                seen_pairs.add(pair_key)
            else:
                print(f"Duplicate exact pair removed in {file}: {mut_name}")

print(len(unique_seqs))

# Save results
output_path = os.path.join(input_dir, f'unique_sequences_{prot_name}.json')
with open(output_path, 'w') as f:
    json.dump(unique_seqs, f, indent=4)
    