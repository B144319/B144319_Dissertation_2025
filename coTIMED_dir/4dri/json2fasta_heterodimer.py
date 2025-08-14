import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def json_to_fasta(json_file, output_fasta):
    with open(json_file, 'r') as f:
        data = json.load(f)

    records = []

    for mut_name, chains in data.items():
        for chain_id, sequence in chains.items():
            seq_id = f"{mut_name}_{chain_id}"
            record = SeqRecord(Seq(sequence), id=seq_id, description="")
            records.append(record)

    SeqIO.write(records, output_fasta, "fasta")
    print(f"FASTA file written to: {output_fasta}")

# Example usage
json_to_fasta("unique_sequences_4dri.json", "4dri_unique_candidates.fasta")


