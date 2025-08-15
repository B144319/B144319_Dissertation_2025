import pandas as pd
import matplotlib.pyplot as plt

# Load the full prediction data
df = pd.read_csv("gnina_summary_combined_postTIMED.csv")

# Extract protein name prefix from filename
df['protein'] = df['filename'].apply(lambda x: x.split('_')[0].lower())

# Split into 2FDE and 4DRI groups
cnn_2fde = df[df['protein'] == '2fde']['CNNscore']
cnn_4dri = df[df['protein'] == '4dri']['CNNscore']

data_to_plot = [cnn_2fde, cnn_4dri]

# Load summary with average CNNaffinity values per protein
#summary_df = pd.read_csv("gnina_docking_summary.csv")
#summary_df['protein'] = summary_df['protein'].str.lower()

# Extract average CNNaffinity values by group
#dots_2fde = summary_df[summary_df['protein'].str.startswith('2fde')]['CNNaffinity_avg']
#dots_4dri = summary_df[summary_df['protein'].str.startswith('4dri')]['CNNaffinity_avg']

# Create the plot
plt.figure(figsize=(8, 6))

# Boxplot
box = plt.boxplot(
    data_to_plot,
    labels=['2FDE-based', '4DRI-based'],
    patch_artist=True,
    medianprops=dict(color='black'),
    showfliers=False
)

# Custom colors for boxes
colors = ['cornflowerblue', 'mediumseagreen']  # 2FDE = blue, 4DRI = green
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

# Overlay average dots
#plt.scatter([1]*len(dots_2fde), dots_2fde, color='navy', edgecolor='white', zorder=3, label='Pre-engineering CNNAffinity (2FDE)')
#plt.scatter([2]*len(dots_4dri), dots_4dri, color='darkgreen', edgecolor='white', zorder=3, label='Pre-engineering CNNAffinity (4DRI)')

# Labels and title
plt.ylabel('CNN Score', fontsize=12)
plt.xlabel('Scaffold', fontsize=12)
plt.title('CNN Score Distribution of Protein Designs', fontsize=14, weight='bold')
plt.axhline(y=0.75, color='red', linestyle='--', linewidth=1.2, label='Cutoff = 0.75')
# Styling
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()

# Save and show
plt.savefig("cnn_score_box_with_avg_dots.png", dpi=300)
plt.show()