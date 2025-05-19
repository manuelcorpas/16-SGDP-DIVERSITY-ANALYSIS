
import matplotlib.pyplot as plt
import pandas as pd

# Ancestry proportions per country (European, Indigenous American, African, Asian)
ancestry_data = {
    'Peru': [25, 65, 8, 2],
    'Mexico': [50, 45, 4, 1],
    'Guatemala': [35, 60, 4, 1],
    'Brazil': [65, 12, 20, 3],
    'Colombia': [60, 25, 12, 3],
    'Chile': [70, 25, 3, 2],
    'Argentina': [80, 15, 3, 2],
    'Uruguay': [85, 10, 3, 2],
    'Dominican Republic': [60, 8, 30, 2],
    'Puerto Rico': [70, 10, 18, 2],
    'Costa Rica': [60, 30, 8, 2],
    'Cuba': [65, 10, 22, 3],
    'Ecuador': [55, 40, 3, 2],
    'Bolivia': [30, 65, 3, 2],
    'Honduras': [45, 40, 12, 3],
    'Panama': [50, 25, 22, 3],
    'Venezuela': [65, 20, 12, 3],
    'Paraguay': [35, 60, 3, 2],
    'El Salvador': [45, 50, 3, 2],
    'Nicaragua': [50, 45, 4, 1],
}

# Define ancestry labels and colors
ancestry_labels = ['European', 'Indigenous American', 'African', 'Asian']
colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd']

# Create DataFrame
df = pd.DataFrame.from_dict(ancestry_data, orient='index', columns=ancestry_labels)

# Optional: sort countries by region or custom order
country_order = [
    'Mexico', 'Guatemala', 'El Salvador', 'Honduras', 'Nicaragua', 'Costa Rica', 'Panama',
    'Cuba', 'Dominican Republic', 'Puerto Rico',
    'Colombia', 'Venezuela', 'Ecuador', 'Peru', 'Bolivia',
    'Brazil', 'Paraguay', 'Chile', 'Argentina', 'Uruguay'
]
df = df.loc[country_order]

# Plot horizontal stacked bar chart
fig, ax = plt.subplots(figsize=(10, 8))
bottom = pd.Series([0]*len(df), index=df.index)

for i, label in enumerate(ancestry_labels):
    ax.barh(df.index, df[label], left=bottom, label=label, color=colors[i])
    bottom += df[label]

# Formatting
ax.set_xlim(0, 100)
ax.set_xlabel('Ancestry Proportion (%)', fontsize=12)
ax.set_title('Ancestry Composition in Latin American Countries', fontsize=14)
ax.invert_yaxis()
ax.legend(loc='upper right', fontsize=10)
ax.grid(axis='x', linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig("latam_ancestry_stacked_bar.png", dpi=600)
plt.savefig("latam_ancestry_stacked_bar.pdf")
plt.show()

