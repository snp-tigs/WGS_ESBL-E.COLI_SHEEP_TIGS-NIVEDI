"""
==============================================================================
Title: Figure 6a - Circos Chord Diagram of Plasmid Replicons
Description: Generates a Circos plot visualizing the shared and unique plasmid 
             replicons among ESBL E. coli isolates from sheep and environmental 
             sources. 
              Script accompanies the publication: Whole-genome analysis reveals 
              shared extended-spectrum β-lactamase-producing Escherichia coli 
              lineages between migratory sheep and surrounding farm environments
==============================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
from pycirclize import Circos
import os

# ==========================================
# 1. Configuration & File Paths
# ==========================================
# Relative paths for GitHub reproducibility
INPUT_CSV   = "./data_main_figures/Figure_6/Fig_6a.csv"

# Output directories and files
RESULTS_DIR = "./results/Figure_6"
OUTPUT_PDF  = os.path.join(RESULTS_DIR, "Fig_6a_Circos_Chord_Diagram.pdf")

# Create results directory if it doesn't exist
os.makedirs(RESULTS_DIR, exist_ok=True)

# ==========================================
# 2. Define Color Aesthetics
# ==========================================
# Sector and Ribbon Colors (Defined directly by Figure 6a Legend)
COLOR_SHEEP_PURPLE = "#8E44AD"  # Deep Purple for Sheep/Rectal Swab
COLOR_ENV_GREEN    = "#27AE60"  # Strong Green for Environment
COLOR_SHARED_YELLOW= "#F1C40F"  # Yellow for shared replicons

# Ribbon colors (RGBA format for transparency)
RIBBON_SHEEP_PURPLE = (0.55, 0.26, 0.67, 0.45) # Purple with 45% opacity
RIBBON_ENV_GREEN    = (0.15, 0.68, 0.37, 0.45) # Green with 45% opacity

district_col = {
    "Chitradurga": "#ff6d3e",
    "Tumkur": "#252d80"
}

village_col = {
    "Dodderi": "#c8755b", "Ganjigunte": "#e18466", "Pavgada ward no 1": "#636ede",
    "Bharamasagara": "#fa9372", "Lokkamanahalli": "#6f7bf7", "Talakere": "#7d88f7",
    "Maddakkanahalli": "#8b95f8", "Pavgada ward no 2": "#9aa2f9",
    "Nerlegudda": "#a8affa", "Anekere": "#b7bdfb", "Herur": "#c5cafb",
    "Dandikere": "#d3d7fc", "Kilarahalli": "#434ca9", "Bagganandu Kaval": "#fa9372",
    "Doddamalligere": "#2f39a0", "Muddapura": "#fa9d80", "Shivapura": "#8288c6",
    "Thirumalapalya": "#5862c5", "Tumkur ward no-9": "#2a3390",
    "Sondekere": "#fbb39c", "Janapanahalli": "#2f39a0", "Kallankere": "#8288c6",
    "Mallapura": "#fcc9b8", "Hosahalli G": "#abafd9", "Vanivilasapura": "#fdded4"
}

stype_col = {
    "Feed": "#adfda2",
    "Soil": "#beb15b",
    "Rectal swab": "#5c3863",
    "Water": "#4d8eff"
}

# ==========================================
# 3. Load, Clean, and Categorize Data
# ==========================================
print("Loading and categorizing data...")
df = pd.read_csv(INPUT_CSV)
df.columns = df.columns.str.strip()

def map_sample_type(val):
    val_str = str(val).strip()
    if 'rectal swab' in val_str.lower():
        return 'Rectal swab'
    return 'Environment'

df['Sample_Category'] = df['Sample Type'].apply(map_sample_type)

def get_exact_subtype(val):
    val_str = str(val).strip()
    for key in stype_col.keys():
        if key.lower() in val_str.lower():
            return key
    return val_str

df['Exact_Subtype'] = df['Sample Type'].apply(get_exact_subtype)
df['GENE'] = df['Plasmid replicon'].astype(str).str.strip()
df['New Sample ID'] = df['New Sample ID'].astype(str).str.strip()
df['District'] = df['District'].astype(str).str.strip()
df['New Village ID'] = df['New Village ID'].astype(str).str.strip()

# Separate and sort samples for clean visual flow
df_env = df[df['Sample_Category'] == 'Environment'].sort_values(
    by=['Sample_Category', 'District', 'New Village ID', 'New Sample ID'])
df_rectal = df[df['Sample_Category'] == 'Rectal swab'].sort_values(
    by=['Sample_Category', 'District', 'New Village ID', 'New Sample ID'])

df_sorted = pd.concat([df_env, df_rectal])

# Aggregate links between Sample and Replicon
link_df = df_sorted.groupby(
    ['New Sample ID', 'Sample_Category', 'Exact_Subtype', 'District', 'New Village ID', 'GENE']
).size().reset_index(name='Weight')

# Mapping dict for sample metadata
sample_meta = df_sorted[
    ['New Sample ID', 'Sample_Category', 'Exact_Subtype', 'District', 'New Village ID']
].drop_duplicates().set_index('New Sample ID').to_dict('index')

# ==========================================
# 4. Classify Replicons (Unique vs Shared)
# ==========================================
gene_categories = df.groupby('GENE')['Sample_Category'].unique()
env_exclusive_genes = []
rectal_exclusive_genes = []
shared_genes = []

for gene, cats in gene_categories.items():
    cat_set = set(cats)
    if cat_set == {'Environment'}:
        env_exclusive_genes.append(gene)
    elif cat_set == {'Rectal swab'}:
        rectal_exclusive_genes.append(gene)
    else:
        shared_genes.append(gene)

env_exclusive_genes.sort()
rectal_exclusive_genes.sort()
shared_genes.sort()

# ==========================================
# 5. Define Sectors & Node Sizes (Optimized Order)
# ==========================================
node_sizes = {}

# A. Shared Replicons (placed centrally to bridge the layout smoothly)
for gene in shared_genes:
    g_weight = link_df[link_df['GENE'] == gene]['Weight'].sum()
    if g_weight > 0: node_sizes[gene] = g_weight

# B. Environment-Exclusive Replicons & Environment Samples
for gene in env_exclusive_genes:
    g_weight = link_df[link_df['GENE'] == gene]['Weight'].sum()
    if g_weight > 0: node_sizes[gene] = g_weight

for sample in df_env['New Sample ID'].unique():
    s_weight = link_df[link_df['New Sample ID'] == sample]['Weight'].sum()
    if s_weight > 0: node_sizes[f"Sample_{sample}"] = s_weight

# C. Rectal Swab Samples & Rectal-Exclusive Replicons
for sample in df_rectal['New Sample ID'].unique():
    s_weight = link_df[link_df['New Sample ID'] == sample]['Weight'].sum()
    if s_weight > 0: node_sizes[f"Sample_{sample}"] = s_weight

for gene in rectal_exclusive_genes:
    g_weight = link_df[link_df['GENE'] == gene]['Weight'].sum()
    if g_weight > 0: node_sizes[gene] = g_weight

# ==========================================
# 6. Initialize Circos & Format Tracks
# ==========================================
print("Building Circos layout...")
circos = Circos(node_sizes, space=0.1)

for sector in circos.sectors:
    sec_name = sector.name

    if sec_name.startswith('Sample_'):
        # --- Format Sample Nodes (Isolates) ---
        sample_id = sec_name.replace('Sample_', '')
        meta = sample_meta.get(sample_id, {})
        exact_type = meta.get('Exact_Subtype', 'Rectal swab')
        dist = meta.get('District', '')
        vill = meta.get('New Village ID', '')

        # Outer Arc (Layer 1): District color
        t_dist = sector.add_track((95, 100))
        t_dist.axis(fc=district_col.get(dist, "#cccccc"), ec="black", lw=0.4)
        t_dist.text(f"ID:{sample_id}", size=5, r=135, orientation="vertical", ha="left")

        # Middle Arc (Layer 2): Village ID color
        t_vill = sector.add_track((90, 94))
        t_vill.axis(fc=village_col.get(vill, "#aaaaaa"), ec="white", lw=0.2)

        # Inner Arc (Layer 3 - closest to center): Sample Type color
        t_base = sector.add_track((85, 89))
        t_base.axis(fc=stype_col.get(exact_type, "#cccccc"), ec="white", lw=0.2)

    else:
        # --- Format Replicon Nodes (Genes) ---
        track = sector.add_track((90, 100))
        
        # Determine sector color based on exclusivity (per Figure Legend)
        if sec_name in shared_genes:
            sector_color = COLOR_SHARED_YELLOW
        elif sec_name in rectal_exclusive_genes:
            sector_color = COLOR_SHEEP_PURPLE
        elif sec_name in env_exclusive_genes:
            sector_color = COLOR_ENV_GREEN
        else:
            sector_color = "lightgrey"
            
        track.axis(fc=sector_color, ec="black", lw=0.5)
        n_count = int(node_sizes[sec_name])
        track.text(f"{sec_name} (n={n_count})", size=7, r=105, orientation="vertical", ha="left", va="center")

# ==========================================
# 7. Draw Links (Chords)
# ==========================================
print("Drawing ribbons...")
current_positions = {node: 0 for node in node_sizes.keys()}

for _, row in link_df.iterrows():
    source = f"Sample_{row['New Sample ID']}"
    target = row['GENE']
    weight = row['Weight']
    cat = row['Sample_Category']

    s_start = current_positions[source]
    s_end = s_start + weight
    current_positions[source] = s_end

    t_start = current_positions[target]
    t_end = t_start + weight
    current_positions[target] = t_end

    # Set ribbon color strictly based on Sample Category
    if cat == 'Rectal swab':
        ribbon_color = RIBBON_SHEEP_PURPLE
    else:
        ribbon_color = RIBBON_ENV_GREEN

    circos.link((source, s_start, s_end), (target, t_start, t_end), color=ribbon_color)

# ==========================================
# 8. Render and Export
# ==========================================
fig = circos.plotfig()
#plt.title("Plasmid Replicon Distribution", fontsize=14, pad=40)

# Save as high-resolution PDF
plt.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight")
print(f"✅ Circos chord diagram successfully saved to: {OUTPUT_PDF}")
