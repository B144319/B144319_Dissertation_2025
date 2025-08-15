import ampal
import pandas as pd
import requests

prots = ['2fdd','2fde','1b6k','6d0d','6d0e','2i4d','6uwc','1yt9','3o9b','6oxp',
         '6w6s','6w6t','6dh1','6dh4','6dh7','3o99','6dgy','5gpg','1fap','4drh',
         '4dri','4drj','3o9a','4m8x','4m8y','2i4w','2i4x','3fap','4fap','3mxe',
         '3gi6','3o9f','3o9g','6oxq','6ode']

prot_info = []

def get_protein_name(pdb_id):
    """Fetch full protein name using RCSB REST API."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        try:
            title = data.get("struct", {}).get("title", "Unknown")
            return title
        except:
            return "Unknown"
    else:
        return "Unknown"

for prot in prots:
    pdb_filename = prot + '.pdb'
    try:
        my_protein = ampal.load_pdb(pdb_filename)
    except Exception as e:
        print(f"Error loading {pdb_filename}: {e}")
        continue

    full_name = get_protein_name(prot)

    chain_num = [0, 1]
    for num in chain_num:
        try:
            mychain = my_protein[num]
            if mychain is not None:
                mychain_id = mychain.id
                chain_weight = mychain.molecular_weight

                prot_info.append({
                    'PDB_ID': prot,
                    'Full_Protein_Name': full_name,
                    'Protein_Molecular_Weight': my_protein.molecular_weight,
                    'Chain_ID': mychain_id,
                    'Chain_Sequence': mychain.sequence,
                    'Chain_Weight': chain_weight,
                    'Chain_Length': len(mychain.sequence)
                })
        except Exception as e:
            print(f"Error processing chain {num} in {prot}: {e}")

prot_df = pd.DataFrame(prot_info)
prot_df.to_csv("candidate_protein_info.csv", index=False)


