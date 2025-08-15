#Script that performs a RDkit similarity search through a database of ligands. 
#Target molecule specified by SMILES structure. 
#Tanimoto similarity score calculated through Morgan Fingerprints.
#Associated proteins from top 20 ranked ligands output to csv file. 

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
import pandas as pd
import ast 


#Read in data from csv file - note that the file is nested: 
df = pd.read_csv('lig_set240715.csv', header=None)

#Convert rows from strings into a dictionary: 
df_dicts = df[0].apply(ast.literal_eval)

#Create flat dictionary from the nested one for easier parsing: 
lig_df = []

for row in df_dicts: 
     for key, val in row.items():
         lig_df.append({
         'ID' : key, 
         'Name' : val['mol_name'], 
         'CACTVS_Smiles' : val['mol_stru']['mol_smile'].get('CACTVS'),
         'OpenEye_Smiles' : val['mol_stru']['mol_smile'].get('OpenEye OEToolkits'),
  })

#convert to pandas dataframe: 
lig_df = pd.DataFrame(lig_df)

#Create a function that will go through the dataframe row by row and create morgan fingerprints
fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=128) #reducing bit length, loses specificity of functional groups but just want a backbone for now. 

def get_molecule_and_fingerprint(smiles): 
    try:
        mol = Chem.MolFromSmiles(smiles)
        fp = fpgen.GetFingerprint(mol)
        return mol, fp
    except:
        return None, None

#Define which column to use for the Smiles format:
smiles_column = 'CACTVS_Smiles'

#Apply this function to each smile in the smiles column, create new columns with molecules and fingerprints: 
lig_df['rdkit_mol'] = lig_df[smiles_column].apply(lambda x: get_molecule_and_fingerprint(x)[0])
lig_df['fingerprint'] = lig_df[smiles_column].apply(lambda x: get_molecule_and_fingerprint(x)[1])

#Mycolactone SMILES format and fingerprint - need to query this against each row in the dataframe: 
myc = 'C[C@H]1C/C(=C/C[C@@H](OC(=O)CCC[C@@H]1OC(=O)/C=C/C(=C/C(=C/C=C/C(=C/[C@@H]([C@H](C[C@H](C)O)O)O)/C)/C)/C)[C@@H](C)C/C(=C/[C@@H](C)[C@@H](C[C@@H](C)O)O)/C)/C'
myc_mol, myc_fp = get_molecule_and_fingerprint(myc)

#Create a function to calculate Tanimoto and Dice similarity metrics from fingerprints: 
def calculate_similarity(row, query_fp): 
    try:
       return (
       DataStructs.TanimotoSimilarity(query_fp, row['fingerprint']), 
       DataStructs.DiceSimilarity(query_fp, row['fingerprint']), 
       DataStructs.CosineSimilarity(query_fp, row['fingerprint'])
       )
    except:
        return (0.0, 0.0, 0.0) #incase of error.

#Add similarity columns to dataframe, applying the similarity function to each row. 
lig_df['Tanimoto_similarity'] = lig_df.apply(lambda row: calculate_similarity(row, myc_fp)[0], axis=1)
lig_df['Dice_similarity'] = lig_df.apply(lambda row: calculate_similarity(row, myc_fp)[1], axis=1)
lig_df['Cosine_similarity'] = lig_df.apply(lambda row: calculate_similarity(row, myc_fp)[2], axis=1)

#create a new dataframe with the similarity scores

results_df = lig_df[['ID', 'Name', smiles_column, 'fingerprint', 'Tanimoto_similarity', 'Dice_similarity', 'Cosine_similarity']].sort_values(by='Tanimoto_similarity', ascending=False)
top10 = results_df['ID'].iloc[0:20].tolist()
import pandas as pd

# Extract top 20 entries
top20_df = results_df[['ID', 'Tanimoto_similarity']].iloc[0:20]

# Save to CSV
top20_df.to_csv('top20_tanimoto.csv', index=False)

all_tanimoto = results_df[['ID', 'Tanimoto_similarity']]

all_tanimoto.to_csv('tanimoto_scores_full.csv')


#Finding binding partners for the most similar ligands: 

#Read csv file, note that it is not in a nested structure and can be used directly. 

prot_df = pd.read_csv("Prot_lig_dataset.csv")

#create a function which matches the ligands in top20 to the proteins in prot_df: 

filtered = []


def match_ligs_to_prots(top_ligs): 
    for lig in top_ligs:
        matches = prot_df[prot_df['ligand'].str.startswith(lig)]
        filtered.append(matches)
    return filtered

           
match_ligs_to_prots(top10)

#convert into one data frame: 
filtered_df = pd.concat(filtered, ignore_index=True)

print (filtered_df['prot_id'], filtered_df['ligand'])

print(filtered_df)

#Save names of matched proteins into a list: 
prot_candidates = filtered_df['prot_id']

prot_candidates = ','.join(prot_candidates.tolist())

#Create a file containing the information of these candidate proteins. 

output_file = open("candidate_proteins.txt", 'w')
output_file.write(prot_candidates)
output_file.close()


