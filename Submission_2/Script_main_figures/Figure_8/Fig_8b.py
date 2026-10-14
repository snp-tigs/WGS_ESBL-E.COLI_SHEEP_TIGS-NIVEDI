"""
==============================================================================
Title: Figure 8b - Arc Diagram of Average Nucleotide Identity (ANI)
Description: Generates an arc diagram showing the ANI of study isolates based 
             on chromosomal contigs (>10,000 bp). Node sizes represent the 
             number of connections (edges), and arc colors depict percent identity.
              Script accompanies the publication: Whole-genome analysis reveals 
              shared extended-spectrum β-lactamase-producing Escherichia coli 
              lineages between migratory sheep and surrounding farm environments
==============================================================================
"""

import os
import csv
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.lines import Line2D

# ==========================================
# 1. Configuration & File Paths
# ==========================================
# Relative paths for GitHub reproducibility
NODES_CSV = "./data_main_figures/Figure_8/Fig_8b_nodes.csv"
EDGES_CSV = "./data_main_figures/Figure_8/Fig_8b_edges.csv"

# Output directories and files
RESULTS_DIR = "./results/Figure_8"
OUTPUT_PDF  = os.path.join(RESULTS_DIR, "Fig_8b_ANI_Arc_Diagram.pdf")

# Create results directory if it doesn't exist
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- COLORS ---
node_colors = {
    "Rectal swab": "#5c3863",
    "Feed": "#adfda2",
    "Soil": "#beb15b",
    "Water": "#4d8eff",
    "water": "#4d8eff" # Fallback for casing
}
default_color_palette = ["#377eb8", "#984ea3", "#ff7f00"]  # Fallbacks for new types

edge_colors = {
    "70-80": "#d9d9d9",
    "80-90": "#ffc872",
    "90-95": "#c957bc",
    "95-100": "#752092"
}

# ==========================================
# 2. Load and Sort Nodes
# ==========================================
print("Loading node metadata...")
nodes = []
seen_types = set()

with open(NODES_CSV, "r") as f:
    delim = '\t' if '\t' in f.readline() else ','
    f.seek(0)
    reader = csv.DictReader(f, delimiter=delim)
    for row in reader:
        # Standardize formatting
        row["Sample_Type"] = str(row["Sample_Type"]).strip()
        nodes.append(row)
        seen_types.add(row["Sample_Type"])

# Assign fallback colors to unknown Sample Types
color_idx = 0
for stype in seen_types:
    if stype not in node_colors:
        node_colors[stype] = default_color_palette[color_idx % len(default_color_palette)]
        color_idx += 1

# SORTING: Sort alphabetically by Sample_Type, then Numerically by Id for structured clustering
nodes.sort(key=lambda x: (x["Sample_Type"], int(x["Id"])))

# Map sorted nodes to X-axis coordinates
node_x = {str(node["Id"]): i for i, node in enumerate(nodes)}

# Initialize degree tracking dictionary
node_degree = {str(node["Id"]): 0 for node in nodes}

# ==========================================
# 3. Load Edges and Calculate Degrees
# ==========================================
print("Loading edge relationships...")
edges = []
with open(EDGES_CSV, "r") as f:
    delim = '\t' if '\t' in f.readline() else ','
    f.seek(0)
    reader = csv.DictReader(f, delimiter=delim)
    for row in reader:
        src = str(row["Source"]).strip()
        tgt = str(row["Target"]).strip()
        weight = float(row["Weight"])

        if src not in node_x or tgt not in node_x or src == tgt:
            continue

        # Categorize sequence identity (ANI)
        if 70.0 <= weight < 80.0:
            cat = "70-80"
        elif 80.0 <= weight < 90.0:
            cat = "80-90"
        elif 90.0 <= weight < 95.0:
            cat = "90-95"
        elif 95.0 <= weight <= 100.0:
            cat = "95-100"
        else:
            continue

        edges.append({"Source": src, "Target": tgt, "Category": cat})

        # Increment true connection degree for both nodes
        node_degree[src] += 1
        node_degree[tgt] += 1

# ==========================================
# 4. Plotting Setup
# ==========================================
print("Generating Arc Diagram...")
fig, ax = plt.subplots(figsize=(16, 8))

# Find min and max degrees for dynamic scaling of actual plot nodes
degrees_list = list(node_degree.values())
min_deg = min(degrees_list) if degrees_list else 0
max_deg = max(degrees_list) if degrees_list else 150
if max_deg == min_deg:
    max_deg = min_deg + 1

# --- Draw Arcs (Edges) ---
for edge in edges:
    src = edge["Source"]
    tgt = edge["Target"]
    cat = edge["Category"]

    x1 = node_x[src]
    x2 = node_x[tgt]

    center_x = (x1 + x2) / 2.0
    dist = abs(x2 - x1)

    arc = Arc((center_x, 0), width=dist, height=dist * 1,
              angle=0, theta1=0, theta2=180,
              edgecolor=edge_colors[cat],
              linewidth=0.5, alpha=0.4, zorder=0.1)
    ax.add_patch(arc)

# --- Draw Nodes & Labels ---
for node in nodes:
    node_id = str(node["Id"])
    x_pos = node_x[node_id]
    stype = node["Sample_Type"]
    c = node_colors[stype]

    # Scale dynamic size based on actual connection degree (min size = 6, max size = 22)
    deg = node_degree.get(node_id, 0)
    dynamic_size = 6 + 16 * ((deg - min_deg) / (max_deg - min_deg))

    # Draw node
    ax.plot(x_pos, -0.5, marker='o', markersize=dynamic_size, color=c, markeredgecolor='black', zorder=2)

    # Draw ID label below
    ax.text(x_pos, -3, node["Id"], rotation=90, ha='center', va='top', fontsize=9)

# ==========================================
# 5. Formatting & Legends
# ==========================================
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_ticks([])

ax.autoscale_view()
plt.subplots_adjust(bottom=0.2)

# Build Dynamic Legend
legend_elements = []

# 1. Sample Type Colors
for stype in sorted(seen_types):
    # Capitalize the label for the legend
    display_label = stype.capitalize() if stype.lower() == 'water' else stype
    legend_elements.append(Line2D([0], [0], marker='o', color='w', label=display_label,
                                  markerfacecolor=node_colors[stype], markeredgecolor='black', markersize=10))

# 2. Edge Categories (Percent Identity)
legend_elements.append(Line2D([0], [0], color='w', label=' '))
legend_elements.append(Line2D([0], [0], color='#752092', lw=4, label='95-100% Identity'))
legend_elements.append(Line2D([0], [0], color='#c957bc', lw=4, label='90-95% Identity'))
legend_elements.append(Line2D([0], [0], color='#ffc872', lw=4, label='80-90% Identity'))
legend_elements.append(Line2D([0], [0], color='#d9d9d9', lw=4, label='70-80% Identity'))

# 3. Node Size (Connections)
legend_elements.append(Line2D([0], [0], color='w', label=' '))
legend_elements.append(Line2D([0], [0], color='w', label='Node Size (Connections):'))

rep_values = [15, 50, 100, 150]
for val in rep_values:
    # Map these fixed values to marker sizes using the same linear scaling logic
    rep_marker_size = 6 + 16 * max(0.0, min(1.0, (val - min_deg) / (max_deg - min_deg)))
    legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                  label=f'{val} connections',
                                  markerfacecolor='gray', markeredgecolor='black', markersize=rep_marker_size))

ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.01, 1), frameon=True, title="Legend")
plt.title("Chromosomal ANI Arc Diagram (>=10kb Contigs)", fontsize=14, pad=20, fontweight='bold')

# ==========================================
# 6. Export Output
# ==========================================
plt.savefig(OUTPUT_PDF, format="pdf", bbox_inches="tight", dpi=300)
print(f"✅ Arc diagram successfully saved as: {OUTPUT_PDF}")
plt.close()
