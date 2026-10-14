"""
==============================================================================
Title: Figures 5b & 5c - Alpha Diversity of Chromosomal Virulence Genes
Description: Calculates and visualizes the Simpson and Shannon alpha diversity 
             indices for chromosomal virulence genes across Sheep and 
             Environmental isolates.
              Script accompanies the publication: Whole-genome analysis reveals 
              shared extended-spectrum β-lactamase-producing Escherichia coli 
              lineages between migratory sheep and surrounding farm environments
==============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys
import os

# ==========================================
# 1. Configuration & File Paths
# ==========================================
# Relative paths for GitHub reproducibility
MATRIX_PATH = "./data_main_figures/Figure_5/Fig_5b_5c.csv"

# Using the shared metadata file from the Figure 4 directory
META_PATH = "./data_main_figures/Figure_4/Fig_4b_4c_5b_5c_Metadata.csv"

# Output directories and files
RESULTS_DIR = "./results/Figure_5"
OUTPUT_CSV = os.path.join(RESULTS_DIR, "Annotated_Diversity_Matrix_Fig5bc.csv")
OUTPUT_PDF_SIMPSON = os.path.join(RESULTS_DIR, "Fig_5b_Simpson_Diversity.pdf")
OUTPUT_PDF_SHANNON = os.path.join(RESULTS_DIR, "Fig_5c_Shannon_Diversity.pdf")

# Create results directory if it doesn't exist
os.makedirs(RESULTS_DIR, exist_ok=True)

# ==========================================
# 2. Load and Prepare Data
# ==========================================
print("Loading datasets...")
try:
    df_matrix = pd.read_csv(MATRIX_PATH)
    df_meta = pd.read_csv(META_PATH)
except FileNotFoundError as e:
    print(f"CRITICAL ERROR: Could not find input file.\n{e}")
    sys.exit(1)

# Bulletproof Merge Strategy: 
# Extract ONLY the first sequence of digits from both IDs (e.g., "Sample-32_S32" -> "32")
df_matrix['Clean_ID'] = df_matrix['Sample_ID'].astype(str).str.extract(r'(\d+)')[0]
df_meta['Clean_ID'] = df_meta['Sample ID'].astype(str).str.extract(r'(\d+)')[0]

# ==========================================
# 3. Calculate Alpha Diversity
# ==========================================
print("Calculating Simpson and Shannon diversity indices for Virulence Genes...")
# Isolate only the numeric gene columns
gene_cols = [col for col in df_matrix.columns if col not in ['Sample_ID', 'Clean_ID']]
df_genes = df_matrix[gene_cols].fillna(0)

# Calculate relative proportions for diversity indices
row_sums = df_genes.sum(axis=1).replace(0, np.nan)  # Prevent division by zero
proportions = df_genes.div(row_sums, axis=0)

# Shannon Index: H = -sum(p_i * ln(p_i))
shannon = -(proportions * np.log(proportions.replace(0, np.nan))).sum(axis=1, skipna=True)

# Simpson Index: D = 1 - sum(p_i^2)
simpson = 1 - (proportions ** 2).sum(axis=1, skipna=True)

df_matrix['Shannon'] = shannon.fillna(0)
df_matrix['Simpson'] = simpson.fillna(0)

# ==========================================
# 4. Merge Metadata and Format Output
# ==========================================
df_merged = pd.merge(df_matrix, df_meta, on='Clean_ID', how='inner')

# Safety Check
if df_merged.empty:
    print("CRITICAL ERROR: Merge failed. No matching Sample IDs between matrix and metadata.")
    sys.exit(1)

# Rename 'Sample_ID' using the standardized 'New Sample ID' from metadata
df_merged['Original_Sample_ID'] = df_merged['Sample_ID']
df_merged['Sample_ID'] = df_merged['New Sample ID']

# Create the binary grouping based on Sample Type (Sheep vs Environment)
df_merged['Group'] = np.where(
    df_merged['Sample Type'].astype(str).str.strip().str.lower() == 'rectal swab',
    'Sheep_isolates',
    'Environmental_isolates'
)

# Clean up and save the fully annotated dataframe
df_merged = df_merged.drop(columns=['Clean_ID'])
df_merged.to_csv(OUTPUT_CSV, index=False)
print(f"✅ Annotated dataset saved to: {OUTPUT_CSV}")

# ==========================================
# 5. Plotting Function
# ==========================================
def plot_diversity(metric_col, ylabel, output_path):
    """
    Generates a boxplot with overlaid strip plot and Kruskal-Wallis statistics.
    """
    sheep_data = df_merged[df_merged['Group'] == 'Sheep_isolates'][metric_col].dropna().values
    env_data = df_merged[df_merged['Group'] == 'Environmental_isolates'][metric_col].dropna().values

    if len(sheep_data) == 0 or len(env_data) == 0:
        print(f"\n⚠️ WARNING: Cannot plot {metric_col}. Missing data for one or both groups!")
        return

    data_series = [sheep_data, env_data]
    series_names = ['Sheep_isolates', 'Environmental_isolates']
    
    # Keeping the same color palette as Figure 4
    box_colors = ['#EDCBA8', '#B8D1AE']

    plt.figure(figsize=(2.6, 4))
    
    # Draw boxplot and stripplot
    sns.boxplot(data=data_series, palette=box_colors, width=0.6, linewidth=1.2)
    sns.stripplot(data=data_series, color='black', jitter=True, size=3.5, alpha=0.7)

    plt.xlabel('Sample Type', fontweight='bold')
    plt.ylabel(ylabel, fontweight='bold')
    plt.xticks(ticks=range(len(series_names)), labels=['Sheep\nIsolates', 'Environmental\nIsolates'])
    
    # Formatting
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # Kruskal-Wallis H-test
    h_stat, p_val = stats.kruskal(*data_series)

    if p_val <= 0.001:
        significance = '***'
    elif p_val <= 0.01:
        significance = '**'
    elif p_val <= 0.05:
        significance = '*'
    else:
        significance = 'ns'

    # Annotation rendering for P-value bracket
    ax = plt.gca()
    y_max = np.max([np.max(series) for series in data_series if len(series) > 0])
    h = (y_max * 0.02) if y_max > 0 else 0.01
    y = y_max + (y_max * 0.05)

    col = 'k'
    x1, x2 = 0, 1
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, c=col)
    ax.text((x1 + x2) * 0.5, y + h, significance, ha='center', va='bottom', color=col, fontsize=12)

    print(f"\n--- {metric_col} Index Statistics ---")
    print(f"Kruskal-Wallis H-statistic: {h_stat:.4f}")
    print(f"P-value: {p_val:.4f}")
    print(f"Significance: {significance}")

    plt.tight_layout()
    
    # Save the figure as a high-resolution PDF
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved plot to: {output_path}")

# ==========================================
# 6. Execute Plots
# ==========================================
print("\nGenerating Figures...")
plot_diversity('Simpson', 'Simpson Diversity Index\n(Chromosomal Virulence Genes)', OUTPUT_PDF_SIMPSON)
plot_diversity('Shannon', 'Shannon Diversity Index\n(Chromosomal Virulence Genes)', OUTPUT_PDF_SHANNON)

print("\nPipeline completed successfully.")
