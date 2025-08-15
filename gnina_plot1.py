import pandas as pd
import matplotlib.pyplot as plt


df1 = pd.read_csv("gnina_docking_summary_nonHIV1.csv")
df1['protein'] = df1['protein'].str.lower()
highlight_df1 = df1[df1['protein'] == '4dri']
other_df1 = df1[df1['protein'] != '4dri']


df2 = pd.read_csv("gnina_docking_summary_.csv")
df2['protein'] = df2['protein'].str.lower()
highlight_df2 = df2[df2['protein'] == '2fde']
other_df2 = df2[df2['protein'] != '2fde']

df1_sorted = pd.concat([other_df1, highlight_df1]).sort_values('CNNaffinity')
df2_sorted = pd.concat([other_df2, highlight_df2]).sort_values('CNNaffinity')


protein_pos_1 = {protein: i for i, protein in enumerate(df1_sorted['protein'])}
protein_pos_2 = {protein: i for i, protein in enumerate(df2_sorted['protein'])}


fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300, sharey=True)

# FKBP proteins
axes[0].scatter([protein_pos_1[p] for p in df1_sorted['protein']],
                df1_sorted['CNNaffinity_avg'],
                color=['red' if p == '4dri' else 'mediumseagreen' for p in df1_sorted['protein']],
                edgecolor='black', s=40, marker='s')

axes[0].set_xticks(range(len(protein_pos_1)))
axes[0].set_xticklabels(df1_sorted['protein'].str.upper(), rotation=90)
axes[0].set_xlabel('Protein', fontsize=11)
axes[0].set_ylabel('CNN Affinity (avg)', fontsize=11)
axes[0].set_title('A. FKBP Family', fontsize=12, weight='bold')
axes[0].grid(axis='y', linestyle='--', alpha=0.5)

# HIV-1 Proteins
axes[1].scatter([protein_pos_2[p] for p in df2_sorted['protein']],
                df2_sorted['CNNaffinity_avg'],
                color=['red' if p == '2fde' else 'cornflowerblue' for p in df2_sorted['protein']],
                edgecolor='black', s=40, marker='s')

axes[1].set_xticks(range(len(protein_pos_2)))
axes[1].set_xticklabels(df2_sorted['protein'].str.upper(), rotation=90)
axes[1].tick_params(axis='y', labelleft=True)
axes[1].set_xlabel('Protein', fontsize=11)
axes[1].set_ylabel('CNN Affinity (avg)', fontsize=11)
axes[1].set_title('B. HIV-1 Protease Family', fontsize=12, weight='bold')
axes[1].grid(axis='y', linestyle='--', alpha=0.5)


plt.tight_layout()
plt.savefig("gnina_scatterplot_protein_vs_affinity.png", dpi=300, bbox_inches='tight')
plt.show()
