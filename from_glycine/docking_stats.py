#Script that outputs an average of the metrics output from GNINA docking. 
#Input is a GNINA summary csv file containing metrics from all docking models. 
import pandas as pd

# Load the GNINA summary file
df = pd.read_csv("gnina_summary.csv")

df['protein'] = df['filename'].apply(lambda x: x.split("_")[0])  # adjust as needed

# Group by protein
grouped = df.groupby("protein")

# Calculate statistics
summary = grouped.agg({
    "minimizedAffinity": ["mean"],
    "CNNscore": ["mean"],
    "CNNaffinity": ["mean"],
    "filename": "count"
}).reset_index()

summary.columns = [
    "protein",
    "minimizedAffinity_avg", 
    "CNNscore_avg", 
    "CNNaffinity_avg", 
    "n_models"
]

# Save
summary.to_csv("gnina_protein_binding_summary.csv", index=False)
print(summary.head())