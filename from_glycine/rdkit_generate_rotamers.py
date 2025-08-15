from rdkit import Chem
from rdkit.Chem import AllChem
import os 

#Creating a new directory for the output sdf files: 
output_dir = "" 
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#Mycolactone SMILES: 
myc = 'C[C@H]1C/C(=C/C[C@@H](OC(=O)CCC[C@@H]1OC(=O)/C=C/C(=C/C(=C/C=C/C(=C/[C@@H]([C@H](C[C@H](C)O)O)O)/C)/C)/C)[C@@H](C)C/C(=C/[C@@H](C)[C@@H](C[C@@H](C)O)O)/C)/C'
mol = Chem.MolFromSmiles(myc)
mol = Chem.AddHs(mol) #Add hydrogens

#Generating 10 conformers: 
AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42, numThreads=0)

#Calculating energies in order to select the top 5 conformers: 
conformer_energies = []
for conf_id in range(mol.GetNumConformers()):
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
    ff.Minimize(maxIts=200)  # Optimize with UFF
    energy = ff.CalcEnergy()
    conformer_energies.append((conf_id, energy))

# Sort by energy and select top 5:
conformer_energies.sort(key=lambda x: x[1])
selected_conf_ids = [x[0] for x in conformer_energies[:5]]

# Save as SDF files:
for i, conf_id in enumerate(selected_conf_ids):
    writer = Chem.SDWriter(os.path.join(output_dir, f"mycolactone_rotamer_{i+1}.sdf"))
    mol.SetProp("_Name", f"Mycolactone_Rotamer_{i+1}")
    writer.write(mol, confId=conf_id)
    writer.close()