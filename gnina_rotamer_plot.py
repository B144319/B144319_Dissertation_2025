import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re


sns.set(style="whitegrid")
df = pd.read_csv("gnina_summary.csv")


def extract_protein(filename):
    return filename.lower().split('_')[0]


def extract_rotamer_id(filename):
    match = re.search(r"rotamer_(\d+)", filename.lower())
    return int(match.group(1)) if match else None

df['protein'] = df['filename'].apply(extract_protein)
df['rotamer'] = df['filename'].apply(extract_rotamer_id)


def plot_bar(ax, protein_name, bar_color):
    data = df[(df['protein'] == protein_name.lower()) & (df['model_id'] == 0)]

    if data.empty:
        ax.set_title(f"No data for {protein_name.upper()}", fontsize=14, weight='bold')
        return

    # Group by rotamer and average CNN affinity
    grouped = data.groupby('rotamer')['CNNaffinity'].mean().sort_index()

    bars = ax.bar(grouped.index, grouped.values, color=bar_color, edgecolor='black', linewidth=0.8)

    max_height = grouped.max()
    ax.axhline(y=max_height, color='firebrick', linestyle='--', linewidth=1.2)

    # Annotate bar values
    for rot, val in grouped.items():
        ax.text(rot, val + 0.03, f"{val:.2f}", ha='center', va='bottom', fontsize=9, fontweight='semibold')

    ax.set_xlabel('Rotamer Number', fontsize=14, fontweight='bold')
    ax.set_ylabel('CNN Affinity', fontsize=14, fontweight='bold')
    ax.set_title(f'{protein_name.upper()} Docked with Mycolactone Rotamers', fontsize=14, weight='bold')
    ax.set_xticks(grouped.index)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.grid(axis='y', linestyle='--', alpha=0.6)
    


fig, axs = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

# 4DRI on the left
plot_bar(axs[0], '4dri', bar_color='mediumseagreen')

# 2FDE on the right
plot_bar(axs[1], '2fde', bar_color='cornflowerblue')

plt.tight_layout()
plt.savefig("rotamer_barchart_4dri_2fde.png", dpi=600, bbox_inches='tight')
plt.show()





