import os
import matplotlib.pyplot as plt
import numpy as np
import json

# Define colors for the bars
iBlue = "#F6F5F1FF"   # intersection
iRed = "#D2848DFF"    # SRS-only
iGreenDark = "#0F85A0FF"  # Hapdiff-only

# Load data from JSON files
with open("sveval/hapdiff_svevals.json") as infile:
    hapdiff_stats = json.load(infile)

# Define the samples
samples = ["M11AO", "M11BG", "M11AJ", "M11AU", "M11AV", "M11AR", "M11BN", "M11B1", "M11B4", "M11BA", "M11BD", 
           "M11BJ", "M11BM", "M11BO", "M11BR", "M11BU", "M11AH", "M11B7", "M11C1", "M11CV", "M11AI"]

# Define the mapping from samples to family IDs
sample_to_family = {
    "M11AO": "RGP_696_3", "M11AU": "RGP_1608_3", "M11AR": "RGP_1219_3", "M11AV": "RGP_25_3", "M11BN": "RGP_1219_4",
    "M11B1": "RGP_1081_3", "M11B4": "RGP_686_3", "M11BA": "RGP_607_3", "M11BD": "RGP_123_3", "M11BG": "RGP_1360_3",
    "M11BJ": "RGP_1238_3", "M11BM": "RGP_1395_3", "M11BO": "RGP_1630_3", "M11BR": "RGP_858_3", "M11BU": "RGP_731_3",
    "M11AH": "RGP_876_3", "M11B7": "RGP_12_3", "M11C1": "RGP_2040_3", "M11CV": "RGP_558_3", "M11AI": "RGP_1770_3",
    "M11AJ": "RGP_1316_3"
}

# Map samples to family IDs for plotting
family_ids = [sample_to_family[sample] for sample in samples]

# Plotting
fig, ax = plt.subplots(figsize=(15, 6))

bar_width = 0.7  # Adjust the bar width as needed
index = np.arange(len(samples))

# Plot bars for Hapdiff
hapdiff_intersections = [hapdiff_stats[sample]['Intersection'] for sample in samples]
hapdiff_srs_only = [hapdiff_stats[sample]['SRS-only'] for sample in samples]
hapdiff_hapdiff_only = [hapdiff_stats[sample]['Hapdiff-only'] for sample in samples]

bar1 = ax.bar(index - bar_width/8, hapdiff_intersections, bar_width, label='Intersection', color=iBlue, edgecolor='black', linewidth=0.2)
bar2 = ax.bar(index - bar_width/8, hapdiff_srs_only, bar_width, bottom=hapdiff_intersections, label='SRS-only', color=iRed, edgecolor='black', linewidth=0.2)
bar3 = ax.bar(index - bar_width/8, hapdiff_hapdiff_only, bar_width, bottom=np.array(hapdiff_intersections) + np.array(hapdiff_srs_only), label='LRS-only', color=iGreenDark, edgecolor='black', linewidth=0.2)

# Write text inside each stacked bar for Hapdiff
for i, (rect1, rect2, rect3) in enumerate(zip(bar1, bar2, bar3)):
    height1 = rect1.get_height()
    height2 = rect2.get_height()
    height3 = rect3.get_height()
    total_height = height1 + height2 + height3
    ax.text(rect1.get_x() + rect1.get_width() / 2, height1 / 2, int(height1), ha='center', va='center', color='black', rotation=90, fontsize=14)
    ax.text(rect2.get_x() + rect2.get_width() / 2, height1 + height2 / 2, int(height2), ha='center', va='center', color='black', rotation=90, fontsize=14)
    ax.text(rect3.get_x() + rect3.get_width() / 2, height1 + height2 + height3 / 2, int(height3), ha='center', va='center', color='black', rotation=90, fontsize=14)

# Add labels, title, and legend
ax.set_xlabel('Probands', fontdict={'fontsize': 14})
ax.set_ylabel('Counts', fontdict={'fontsize': 14})
ax.set_title('SV Comparison between LRS and SRS', fontdict={'fontsize': 20})
ax.set_xticks(index)
ax.tick_params(axis='y', labelsize=13)
ax.set_ylim(0, 40000)
ax.set_xticklabels(samples, rotation=30, ha='right', fontdict={'fontsize': 13})
ax.legend(loc='upper center', ncol=4, fontsize=14)

plt.savefig('Figure-5B.pdf', dpi=600, format='pdf')
plt.show()
