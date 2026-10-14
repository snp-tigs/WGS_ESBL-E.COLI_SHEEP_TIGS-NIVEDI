"""
==============================================================================
Title: Figure 5e - 3D PCoA of Virulence Genes
Description: Performs Classical Principal Coordinates Analysis (PCoA) using 
             Jaccard distances to visualize the structural divergence of 
             virulence gene profiles between chromosomal and plasmid contigs, 
             stratified by environmental and host (sheep) origin.
             Script accompanies the publication: Whole-genome analysis reveals 
              shared extended-spectrum β-lactamase-producing Escherichia coli 
              lineages between migratory sheep and surrounding farm environments
==============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
import os

# ==========================================
# 1. Configuration & File Paths
# ==========================================
# Relative paths for GitHub reproducibility
CHR_FILE     = "./data_main_figures/Figure_5/Fig_5e_Chromosome_Vir_gene.csv"
PLASMID_FILE = "./data_main_figures/Figure_5/Fig_5e_Plasmid_Vir_gene.csv"

# Output directories and files
RESULTS_DIR      = "./results/Figure_5"
INTERMEDIATE_CSV = os.path.join(RESULTS_DIR, "Fig_5e_Intermediate_PresenceAbsence_Matrix.csv")
OUTPUT_PLOT      = os.path.join(RESULTS_DIR, "Fig_5e_3D_PCoA.pdf") # Saved as high-res PDF for publication

# Create results directory if it doesn't exist
os.makedirs(RESULTS_DIR, exist_ok=True)

# ==========================================
# 2. Load and Prepare Data
# ==========================================
print("Loading and preparing virulence gene data...")
df_chr = pd.read_csv(CHR_FILE)
df_plasmid = pd.read_csv(PLASMID_FILE)

# Strip whitespace from columns
df_chr.columns = df_chr.columns.str.strip()
df_plasmid.columns = df_plasmid.columns.str.strip()

# Standardize gene column names (Ensuring 'Gene' or 'Product of Gene' are mapped correctly)
if 'Product of Gene' in df_plasmid.columns:
    df_plasmid = df_plasmid.rename(columns={'Product of Gene': 'Gene'})
if 'Product of Gene' in df_chr.columns:
    df_chr = df_chr.rename(columns={'Product of Gene': 'Gene'})

# Define the genomic location for each dataset
df_chr['Genomic_Location'] = 'Chromosome'
df_plasmid['Genomic_Location'] = 'Plasmid'

# Extract essential columns and combine
df_chr_sub = df_chr[['Sample ID', 'Sample Type', 'Genomic_Location', 'Gene']]
df_plasmid_sub = df_plasmid[['Sample ID', 'Sample Type', 'Genomic_Location', 'Gene']]

df_combined = pd.concat([df_chr_sub, df_plasmid_sub], ignore_index=True)

# Clean string values
for col in ['Sample ID', 'Sample Type', 'Gene']:
    df_combined[col] = df_combined[col].astype(str).str.strip()

# ==========================================
# 3. Create Binary Presence/Absence Matrix
# ==========================================
print("Generating binary presence/absence matrix...")
# Remove duplicates to ensure strict presence/absence
df_unique = df_combined.drop_duplicates()

# Create composite row IDs
row_ids = df_unique['Sample ID'] + "||" + df_unique['Sample Type'] + "||" + df_unique['Genomic_Location']
df_unique = df_unique.copy()
df_unique['Row_ID'] = row_ids

# Pivot into a wide presence/absence matrix (0 or 1)
df_matrix = df_unique.pivot_table(
    index=['Row_ID', 'Sample ID', 'Sample Type', 'Genomic_Location'],
    columns='Gene',
    values='Gene',
    aggfunc=len,
    fill_value=0
)

# Convert counts to strict binary (1 or 0)
df_binary = (df_matrix > 0).astype(int)

# Reset index to separate metadata for plotting
df_meta = df_binary.index.to_frame(index=False)
abundance_data = df_binary.values

# Save intermediate matrix for transparency and downstream use
df_intermediate = pd.concat([df_meta, pd.DataFrame(abundance_data, columns=df_binary.columns)], axis=1)
df_intermediate.to_csv(INTERMEDIATE_CSV, index=False)
print(f"✅ Intermediate matrix saved to: {INTERMEDIATE_CSV}")

# ==========================================
# 4. Calculate Exact Jaccard Distance & PCoA
# ==========================================
print("Calculating Jaccard distance matrix...")
# pdist with 'jaccard' computes the exact binary Jaccard dissimilarity coefficient
jaccard_condensed = pdist(abundance_data, metric='jaccard')
jaccard_dist = squareform(jaccard_condensed)

print("Performing Classical Principal Coordinates Analysis (PCoA)...")
# Classical MDS (PCoA) matches R's cmdscale()
pcoa = MDS(n_components=3, dissimilarity='precomputed', random_state=1, metric=True)
pcoa_results = pcoa.fit_transform(jaccard_dist)

# Calculate exact Eigenvalues and Percentage of Variance Explained per Axis
n = jaccard_dist.shape[0]
H = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * H.dot(jaccard_dist ** 2).dot(H)
eigenvalues = np.linalg.eigvalsh(B)

eigenvalues = np.sort(eigenvalues)[::-1]  # Sort descending
total_variance = np.sum(eigenvalues[eigenvalues > 0])
variance_explained = (eigenvalues[:3] / total_variance) * 100

print("\n📊 Variance Explained by PCoA Axes:")
print(f"   • Axis 1: {variance_explained[0]:.2f}%")
print(f"   • Axis 2: {variance_explained[1]:.2f}%")
print(f"   • Axis 3: {variance_explained[2]:.2f}%")

# Format plotting dataframe
df_plot = df_meta.copy()
df_plot['PCoA1'] = pcoa_results[:, 0]
df_plot['PCoA2'] = pcoa_results[:, 1]
df_plot['PCoA3'] = pcoa_results[:, 2]

# Standardize sample types for the legend
df_plot['Group'] = np.where(
    df_plot['Sample Type'].str.lower() == 'rectal swab',
    'Sheep Isolates',
    'Environmental Isolates'
)

# ==========================================
# 5. 3D Visualization
# ==========================================
print("\nGenerating 3D plot...")
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

colors = {
    'Chromosome': '#4C72B0',  # Muted Blue
    'Plasmid': '#C44E52'      # Muted Red
}

markers = {
    'Sheep Isolates': 'o',          # Circle
    'Environmental Isolates': '^'   # Triangle
}

for loc in ['Chromosome', 'Plasmid']:
    for grp in ['Sheep Isolates', 'Environmental Isolates']:
        subset = df_plot[(df_plot['Genomic_Location'] == loc) & (df_plot['Group'] == grp)]
        if not subset.empty:
            ax.scatter(
                subset['PCoA1'], subset['PCoA2'], subset['PCoA3'],
                c=colors[loc],
                marker=markers[grp],
                s=120,
                edgecolor='black',
                alpha=0.8,
                label=f"{loc} ({grp})"
            )

# Styling with calculated variance percentages on axis labels
ax.set_title('3D PCoA of Virulence Gene Profiles (Jaccard Distance)', fontsize=14, fontweight='bold')
ax.set_xlabel(f'PCoA Axis 1 ({variance_explained[0]:.1f}%)', fontsize=11)
ax.set_ylabel(f'PCoA Axis 2 ({variance_explained[1]:.1f}%)', fontsize=11)
ax.set_zlabel(f'PCoA Axis 3 ({variance_explained[2]:.1f}%)', fontsize=11)

ax.view_init(elev=20, azim=45)

# Style the legend
ax.legend(title='Genomic Location & Origin', 
          bbox_to_anchor=(1.05, 1), 
          loc='upper left',
          frameon=True)

plt.tight_layout()

# Save the plot as a high-resolution PDF for publication
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"✅ 3D PCoA Plot successfully saved to: {OUTPUT_PLOT}")
