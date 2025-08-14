#must be in the cotimed/aposteriori conda env and must be in the directory containing these files.
#command returns a binary dataset to be fed into cotimed based on a pdb file containing protein 
#of choice with the ligand of choice, and a ligand.csv file containing the name of the ligand 
#in the pdb file e.g. MYC. 

make-frame-dataset . \
    -e .pdb \
    -g True \
    --keep_side_chain_portion 0.5 \
    -ae BackCBSideOrg \
    --voxels-per-side 21 \
    --frame-edge-length 21 \
    -n 2fde_MYC_pose3 \
    --cfile 