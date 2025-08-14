import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def json_to_fasta_homodimer(json_file, output_fasta):
    with open(json_file, 'r') as f:
        data = json.load(f)

    records = []
    for mut_name, sequence in data.items():
        # Create chain A
        record_A = SeqRecord(Seq(sequence), id=f"{mut_name}_chain_A", description="")
        # Create chain B (identical)
        record_B = SeqRecord(Seq(sequence), id=f"{mut_name}_chain_B", description="")
        records.extend([record_A, record_B])

    SeqIO.write(records, output_fasta, "fasta")
    print(f"FASTA homodimer file written to: {output_fasta}")

# Example usage
json_to_fasta_homodimer("unique_sequences_2fde.json", "2fde_unique_candidates_homodimer.fasta")