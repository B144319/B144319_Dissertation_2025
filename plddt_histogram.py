
import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
csv_file = "Boltz_scores_full_combined_postIE.csv"
data = pd.read_csv(csv_file)

# Extract pLDDT values
plddt_values = data['complex_plddt']

# global plot style
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 300,
    'axes.linewidth': 1,
})

# histogram
fig, ax = plt.subplots(figsize=(7, 6))

n, bins, patches = ax.hist(
    plddt_values, bins=10, edgecolor='black',
    color='lightslategrey', rwidth=0.95, zorder=2
)

# Add pLDDT threshold line
ax.axvline(x=0.88, color='red', linestyle='--', linewidth=1.8, label='Cutoff: pLDDT = 0.88')


ax.set_title('Structural Confidence (pLDDT) Distribution of Candidates\nAfter Protein-Protein Interface Engineering', pad=12, fontweight='bold')
ax.set_xlabel('Predicted Local Distance Difference Test (pLDDT)')
ax.set_ylabel('Frequency')
ax.set_xlim(0.7, 1.0)
ax.set_ylim(0, max(n) + 1)
ax.grid(axis='y', linestyle='--', alpha=0.4)
ax.legend(loc='upper left', frameon=False)
plt.tight_layout()
plt.savefig('plddt_histogram_postIE.png')
plt.show()
