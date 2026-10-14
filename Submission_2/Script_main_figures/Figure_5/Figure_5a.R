# ==============================================================================
# Title: Figure 5a - Heatmap of Virulence Genes
# Description: Generates a presence/absence heatmap comparing virulence gene 
#              distribution across chromosomal and plasmid contigs of ESBL E. coli 
#              isolates from sheep and environmental reservoirs.
#              Script accompanies the publication: Whole-genome analysis reveals 
#              shared extended-spectrum β-lactamase-producing Escherichia coli 
#              lineages between migratory sheep and surrounding farm environments
# ==============================================================================

# Load required packages silently
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ==========================================
# 1. Configuration & File Paths
# ==========================================
# Relative paths for GitHub reproducibility
INPUT_CSV   <- "./data_main_figures/Figure_5/Fig_5a_data.csv"

# Output directories and files
RESULTS_DIR <- "./results/Figure_5"
OUTPUT_PDF  <- paste0(RESULTS_DIR, "/Fig_5a_Heatmap.pdf")
OUTPUT_TXT  <- paste0(RESULTS_DIR, "/Fig_5a_Clustered_Columns.txt")

# Create results directory if it doesn't exist
if(!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ==========================================
# 2. Load and Prepare Data
# ==========================================
cat("Loading virulence gene data...\n")
df <- read.csv(INPUT_CSV, check.names = FALSE)

# Sort Data: Group by Sample Type first, then by District
temp_group <- factor(
  ifelse(df$`Sample Type` == "Rectal swab", "Rectal swab", "Environment"), 
  levels = c("Rectal swab", "Environment")
)
df <- df[order(temp_group, df$District), ]

# Extract metadata (Columns 1-4) and genomic matrix (Columns 5 onwards)
metadata <- df[, 1:4]
vir_matrix <- as.matrix(df[, 5:ncol(df)])

# Set row names to Sample ID to label the heatmap y-axis
rownames(vir_matrix) <- metadata$`Sample ID`

# ==========================================
# 3. Define Visual Aesthetics & Palettes
# ==========================================
district_col <- c(
  "Chitradurga" = "#ff6d3e",
  "Tumkur"      = "#252d80"
)

village_col <- c(
  "Dodderi"           = "#c8755b", "Ganjigunte"        = "#e18466",
  "Pavgada ward no 1" = "#636ede", "Bharamasagara"     = "#fa9372",
  "Lokkamanahalli"    = "#6f7bf7", "Talakere"          = "#7d88f7",
  "Maddakkanahalli"   = "#8b95f8", "Pavgada ward no 2" = "#9aa2f9",
  "Nerlegudda"        = "#a8affa", "Anekere"           = "#b7bdfb",
  "Herur"             = "#c5cafb", "Dandikere"         = "#d3d7fc",
  "Kilarahalli"       = "#434ca9", "Bagganandu Kaval"  = "#fa9372",
  "Doddamalligere"    = "#2f39a0", "Muddapura"         = "#fa9d80",
  "Shivapura"         = "#8288c6", "Thirumalapalya"    = "#5862c5",
  "Tumkur ward no-9"  = "#2a3390", "Sondekere"         = "#fbb39c",
  "Janapanahalli"     = "#2f39a0", "Kallankere"        = "#8288c6",
  "Mallapura"         = "#fcc9b8", "Hosahalli G"       = "#abafd9",
  "Vanivilasapura"    = "#fdded4"
)

stype_col <- c(
  "Feed"        = "#adfda2",
  "Soil"        = "#beb15b",
  "Rectal swab" = "#5c3863",
  "Water"       = "#4d8eff"
)

# Presence/absence categorical colors strictly matched to the figure legend
presence_colors <- c(
  "0" = "#E0E0E0",   # Grey: Absent
  "1" = "#00B1D2",   # Light Blue: Present in chromosome
  "2" = "#EF8A62",   # Peach: Present in plasmids
  "3" = "#343148"    # Dark Blue: Present in both
)

# ==========================================
# 4. Create Heatmap Annotations
# ==========================================
# Grouping variable for structural row splits
row_group <- factor(
  ifelse(metadata$`Sample Type` == "Rectal swab", "Rectal swab", "Environment"), 
  levels = c("Rectal swab", "Environment")
)

# Left-side metadata annotation bars
row_annot <- rowAnnotation(
  District   = metadata$District,
  Village    = metadata$Village, SampleType = metadata$`Sample Type`,
  col = list(
    District   = district_col,
    Village    = village_col,
    SampleType = stype_col
  ),
  annotation_name_side = "top",
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    District   = list(title = "District"),
    Village    = list(title = "Village"),
    SampleType = list(title = "Sample Type")
  )
)

# Custom legend for the main matrix
custom_legend <- Legend(
  labels = c("Absent", "Present in chromosome", "Present in plasmids", "Present in both"),
  title = "Gene Location",
  legend_gp = gpar(fill = c("#E0E0E0", "#00B1D2", "#EF8A62", "#343148"))
)

# ==========================================
# 5. Build and Export Heatmap
# ==========================================
ht <- Heatmap(
  vir_matrix,
  show_heatmap_legend = FALSE, 
  col = presence_colors,
  
  # Disable hierarchical clustering to strictly respect manual Sample ID sorting
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # Apply the custom row split separating Rectal swab and Environment
  row_split = row_group,
  row_title = "%s", 
  
  show_row_names = TRUE, 
  show_column_names = TRUE,
  column_names_side = "bottom",
  row_names_side = "left",
  left_annotation = row_annot,
  rect_gp = gpar(col = "black", lwd = 0.5)
)

cat("Generating and saving high-resolution PDF...\n")
pdf(OUTPUT_PDF, width = 16, height = 10)
ht_drawn <- draw(ht,
                 annotation_legend_list = list(custom_legend),
                 heatmap_legend_side = "right",
                 annotation_legend_side = "right"
)
dev.off()

# ==========================================
# 6. Extract Clustered Column Order
# ==========================================
col_indices <- column_order(ht_drawn)

if(is.list(col_indices)) {
  col_indices <- unlist(col_indices)
}

clustered_col_names <- colnames(vir_matrix)[col_indices]

# Save column order to text file for downstream use
writeLines(clustered_col_names, OUTPUT_TXT)

cat("✅ Success! Files saved to:\n")
cat(" - PDF:", OUTPUT_PDF, "\n")
cat(" - TXT:", OUTPUT_TXT, "\n")
