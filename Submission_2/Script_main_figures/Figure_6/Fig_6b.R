# ==============================================================================
# Title: Figure 6b - Plasmid Replicon & Gene Co-occurrence Heatmap
# Description: Generates a presence/absence heatmap mapping the co-occurrence of 
#              plasmid replicons with both virulence genes and ARGs. Features
#              row/column hierarchical clustering with stacked payload bar plots.
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
})

# ==========================================
# 1. Configuration & Data Loading
# ==========================================
# Relative paths for GitHub reproducibility
INPUT_CSV   <- "./data_main_figures/Figure_6/Fig_6b.csv"

# Output directories and files
RESULTS_DIR <- "./results/Figure_6"
OUTPUT_PDF  <- paste0(RESULTS_DIR, "/Fig_6b_Cooccurrence_Heatmap.pdf")
OUTPUT_CSV  <- paste0(RESULTS_DIR, "/Fig_6b_Heatmap_Gene_Order.csv")

# Create results directory if it doesn't exist
if(!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# Read the data (safely handle row names or ID column mapping)
cat("Loading data...\n")
df <- read.csv(INPUT_CSV, check.names = FALSE, stringsAsFactors = FALSE)

if ("Plasmid" %in% colnames(df)) {
  rownames(df) <- df$Plasmid
  df$Plasmid <- NULL
} else {
  # Safe Fallback: Assume first column is Plasmid ID
  rownames(df) <- df[, 1]
  df <- df[, -1]
}

# Convert to strict numeric binary matrix
mat <- as.matrix(df)
mat[is.na(mat)] <- 0
mat <- ifelse(mat > 0, 1, 0)

# ==========================================
# 2. Gene Classification (ARG vs Virulence)
# ==========================================
# Exact AMR gene checklist for classification
args_list <- c("CmlA7", "aac(3)-IId", "aac(6')-Ib-cr", "aadA2", "aadA5", 
               "ant(3'')-Ia", "aph(3'')-Ib", "aph(3')-Ia", "aph(6)-Id", 
               "blaCMY-2", "blaCTX-M-15", "blaCTX-M-55", "blaDHA-1", 
               "blaOXA-1", "blaSHV-12", "blaTEM-104", "blaTEM-1B", 
               "blaTEM-30", "dfrA1", "dfrA12", "dfrA14", "dfrA17", 
               "dfrA5", "dfrA7", "mef(B)", "mph(A)", "qnrB4", "qnrB7", 
               "qnrS1", "qnrS10", "qnrS13", "sul1", "sul2", "sul3", "tet(A)")

# Create metadata factor: Classify each column based on search parameters
gene_types <- ifelse(colnames(mat) %in% args_list, "ARG", "Virulence")

# Standard color palettes (Strictly mapped to Figure 6b legend requirements)
type_colors <- c("ARG" = "#6A6294", "Virulence" = "#FBBD41") # Purple & Yellow

# ==========================================
# 3. Calculate Frequencies for Bar Charts
# ==========================================
# Total replicon count for each individual gene (X-axis/horizontal)
gene_frequencies <- colSums(mat)

# Mathematical matrix splitting to sum loads by category (Y-axis/vertical)
arg_counts <- rowSums(mat[, gene_types == "ARG", drop = FALSE])
vir_counts <- rowSums(mat[, gene_types == "Virulence", drop = FALSE])

# Bind loads together for row-stacked annotation
plasmid_payloads <- cbind(ARG = arg_counts, Virulence = vir_counts)

# ==========================================
# 4. Define Annotations & Custom Legends
# ==========================================
# Top Annotation: Gene counts and Group colored strip (ARG vs Virulence)
top_anno <- HeatmapAnnotation(
  Gene_Frequency = anno_barplot(
    gene_frequencies, 
    gp = gpar(fill = "#525252", border = NA), # Dark grey bars
    height = unit(2, "cm"),
    axis_param = list(side = "left", labels_rot = 0)
  ),
  Gene_Type = gene_types,
  col = list(Gene_Type = type_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

# Right Annotation: Stacked bar chart tracing plasmid load
right_anno <- rowAnnotation(
  Gene_Load = anno_barplot(
    plasmid_payloads, 
    gp = gpar(fill = c("#6A6294", "#FBBD41"), border = NA), # Purple & Yellow
    width = unit(3, "cm"),
    axis_param = list(side = "bottom")
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

# Legend definitions mapping annotations to graphic label specs
payload_legend <- Legend(
  labels = c("ARG Count", "Virulence Count"),
  legend_gp = gpar(fill = c("#6A6294", "#FBBD41")),
  title = "Plasmid Payload"
)

# Matrix presence-absence colors
heatmap_cols <- c("0" = "#f5f5f5", "1" = "#1f78b4") # Grey & Light Blue

# ==========================================
# 5. Build, Clustering, and Export Heatmap
# ==========================================
ht <- Heatmap(
  mat,
  name = "Presence",
  col = heatmap_cols,
  column_split = gene_types,             # Separate visual clusters
  cluster_rows = TRUE,                   # Euclidean H-clustering
  cluster_columns = TRUE, 
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),   
  column_names_gp = gpar(fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 0.5), # Grid boundaries
  top_annotation = top_anno,
  right_annotation = right_anno,
  row_title = "Plasmid Replicons",
  column_title = "Gene Type",
  heatmap_legend_param = list(
    title = "Presence",
    labels = c("Absent", "Present"),
    at = c(0, 1)
  )
)

cat("Drawing and saving heatmap PDF...\n")
pdf(OUTPUT_PDF, width = 16, height = 10)
ht_drawn <- draw(
  ht, 
  merge_legend = TRUE, 
  heatmap_legend_side = "right", 
  annotation_legend_side = "right",
  annotation_legend_list = list(payload_legend)
)
dev.off()

# ==========================================
# 6. Extract Clustered Order
# ==========================================
# Extract column order index arrays from the finalized plot
col_order_indices <- unlist(column_order(ht_drawn))
ordered_genes <- colnames(mat)[col_order_indices]

# Save left-to-right ordered genes to CSV
write.csv(data.frame(Clustered_Gene_Order = ordered_genes), 
          OUTPUT_CSV, 
          row.names = FALSE)

cat("✅ Success! Files saved to:\n")
cat(" - PDF:", OUTPUT_PDF, "\n")
cat(" - CSV:", OUTPUT_CSV, "\n")
