import matplotlib.pyplot as plt
import numpy as np

# Data
metrics = ["CNN Affinity", "Minimized Affinity\n(Vina Score, kcal/mol)", "CNN Score"]
original_values = [8.35, -9.69, 0.58]
final_values = [8.62, -10.49, 0.81]

x = np.arange(len(metrics))  
width = 0.35  

# Create plot
fig, ax = plt.subplots(figsize=(7, 5), dpi=300)

bars1 = ax.bar(x - width/2, original_values, width, label='4DRI (Scaffold)', color="mediumseagreen")
bars2 = ax.bar(x + width/2, final_values, width, label='4dri_p324c1 (Final Design)', color="slategrey")

# Labels above bars
def label_bars(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

label_bars(bars1)
label_bars(bars2)


ax.set_ylabel('Value', fontsize=11)
ax.set_title('GNINA Docking Metrics', fontsize=13, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(metrics, fontsize=10)
ax.legend(frameon=False, fontsize=9)
ax.axhline(0, color='black', linewidth=0.8)  
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.tight_layout()


plt.savefig('gnina_final_bar.png')
plt.show()
