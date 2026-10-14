# ==============================================================================
# Title: Figure 4d - Chao 1 Alpha Diversity of ARGs
# Description: Calculates and visualizes the non-parametric Chao 1 alpha diversity 
#              index for antimicrobial resistance genes (ARGs) across global 
#              chromosomal/plasmid contigs and stratified sheep/environmental ESBL E.coli isolates.
#              Script accompanies the publication: Whole-genome analysis reveals 
#              shared extended-spectrum β-lactamase-producing Escherichia coli 
#              lineages between migratory sheep and surrounding farm environments
# ==============================================================================

# Load required packages silently
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
})

# ==========================================
# 1. Configuration & File Paths
# ==========================================
# Relative paths for GitHub reproducibility
CHR_CSV     <- "./data_main_figures/Figure_4/Fig_4d_Chromosome_ARGs.csv"
PLASMID_CSV <- "./data_main_figures/Figure_4/Fig_4d_Plasmid_ARGs.csv"

# Output directories and files
RESULTS_DIR <- "./results/Figure_4"
OUTPUT_PLOT <- paste0(RESULTS_DIR, "/Fig_4d_Chao1_Alpha_Diversity.pdf")

# Create results directory if it doesn't exist
if(!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ==========================================
# 2. Load and Prepare Data
# ==========================================
message("Loading chromosomal and plasmid ARG datasets...")
df_chr <- read.csv(CHR_CSV, check.names = FALSE, stringsAsFactors = FALSE)
df_plasmid <- read.csv(PLASMID_CSV, check.names = FALSE, stringsAsFactors = FALSE)

# Clean column names and Sample IDs
colnames(df_chr) <- trimws(colnames(df_chr))
colnames(df_plasmid) <- trimws(colnames(df_plasmid))
df_chr$`Sample ID` <- trimws(as.character(df_chr$`Sample ID`))
df_plasmid$`Sample ID` <- trimws(as.character(df_plasmid$`Sample ID`))

meta_cols <- c("Sample ID", "District", "Village", "Sample Type")

# Extract gene matrices
chr_gene_cols <- setdiff(colnames(df_chr), meta_cols)
plasmid_gene_cols <- setdiff(colnames(df_plasmid), meta_cols)

mat_chr <- as.matrix(df_chr[, chr_gene_cols, drop = FALSE])
mat_plasmid <- as.matrix(df_plasmid[, plasmid_gene_cols, drop = FALSE])

# Ensure strict binary values (0 or 1)
mat_chr[is.na(mat_chr)] <- 0
mat_chr <- ifelse(mat_chr > 0, 1, 0)

mat_plasmid[is.na(mat_plasmid)] <- 0
mat_plasmid <- ifelse(mat_plasmid > 0, 1, 0)

# ==========================================
# 3. Calculate Alpha Diversity & Formats
# ==========================================
# Calculate Chao 1 Diversity (equivalent to row sums for presence/absence)
df_chr$Chao_Chr <- rowSums(mat_chr)
df_plasmid$Chao_Plasmid <- rowSums(mat_plasmid)

# Merge diversity metrics by Sample ID
df_merged <- merge(
  df_chr[, c("Sample ID", "Sample Type", "Chao_Chr")],
  df_plasmid[, c("Sample ID", "Chao_Plasmid")],
  by = "Sample ID",
  all = TRUE
)
df_merged[is.na(df_merged)] <- 0

# Standardize Environmental vs Host Sources
standardize_source <- function(stype) {
  stype_clean <- trimws(tolower(stype))
  if (stype_clean == "rectal swab") {
    return("Rectal Swab")
  } else if (stype_clean %in% c("soil", "feed", "water")) {
    return("Environment")
  } else {
    return(NA)
  }
}
df_merged$Source_Group <- sapply(df_merged$`Sample Type`, standardize_source)

# --- Dataset 1: Global (Chromosomal vs Plasmid) ---
df_global <- df_merged %>%
  pivot_longer(cols = c(Chao_Chr, Chao_Plasmid), 
               names_to = "Compartment", 
               values_to = "Alpha_Diversity") %>%
  mutate(
    Group = recode(Compartment, 
                   "Chao_Chr" = "Chromosomal (Global)", 
                   "Chao_Plasmid" = "Plasmid (Global)"),
    Point_Color = "black", # Explicitly forced to black for global boxes
    Point_Shape = 16       # 16 = Solid Circle
  ) %>%
  select(`Sample ID`, Group, Alpha_Diversity, Point_Color, Point_Shape)

# --- Dataset 2: Stratified (4-Group Breakdown) ---
df_chr_strat <- df_merged %>%
  filter(!is.na(Source_Group)) %>%
  select(`Sample ID`, Source_Group, Chao_Chr) %>%
  mutate(
    Group = paste0("Chr - ", Source_Group), 
    Alpha_Diversity = Chao_Chr, 
    Point_Color = "#4C72B0",
    Point_Shape = ifelse(Source_Group == "Rectal Swab", 16, 17)
  )

df_plas_strat <- df_merged %>%
  filter(!is.na(Source_Group)) %>%
  select(`Sample ID`, Source_Group, Chao_Plasmid) %>%
  mutate(
    Group = paste0("Plasmid - ", Source_Group), 
    Alpha_Diversity = Chao_Plasmid, 
    Point_Color = "#C44E52",
    Point_Shape = ifelse(Source_Group == "Rectal Swab", 16, 17)
  )

df_stratified <- bind_rows(
  df_chr_strat %>% select(`Sample ID`, Group, Alpha_Diversity, Point_Color, Point_Shape),
  df_plas_strat %>% select(`Sample ID`, Group, Alpha_Diversity, Point_Color, Point_Shape)
)

# Combine all 6 categories into a single unified dataframe
df_all_six <- bind_rows(df_global, df_stratified)

# Set explicit factor ordering for all 6 boxes on the single x-axis
df_all_six$Group <- factor(df_all_six$Group, levels = c(
  "Chromosomal (Global)", 
  "Plasmid (Global)",
  "Chr - Rectal Swab", 
  "Chr - Environment",
  "Plasmid - Rectal Swab", 
  "Plasmid - Environment"
))

# Compute maximum y-value for proper bracket headroom
max_y <- max(df_all_six$Alpha_Diversity, na.rm = TRUE)

# ==========================================
# 4. Visualization (Unified 6-Box Plot)
# ==========================================
message("Generating unified 6-box plot...")

# Background box colors
box_colors <- c(
  "Chromosomal (Global)"  = "#00B1D2",
  "Plasmid (Global)"      = "#b86544",
  "Chr - Rectal Swab"     = "#F4A460", 
  "Chr - Environment"     = "#B8D1AE", 
  "Plasmid - Rectal Swab" = "#F4A460", 
  "Plasmid - Environment" = "#8FBC8F"
)

# Define comparison brackets across the boxes
comparisons_list <- list(
  c("Chromosomal (Global)", "Plasmid (Global)"),       # 1. Global comparison
  c("Chr - Rectal Swab", "Chr - Environment"),         # 2. Chr source comparison
  c("Chr - Rectal Swab", "Plasmid - Rectal Swab"),     # 3. Rectal compartment comparison
  c("Plasmid - Rectal Swab", "Plasmid - Environment"), # 4. Plasmid source comparison
  c("Plasmid - Environment", "Chr - Environment")      # 5. Environment compartment comparison
)

p <- ggplot(df_all_six, aes(x = Group, y = Alpha_Diversity, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.8, whisker.linewidth = 0.12, 
               median.linewidth = 0.001, box.linetype = 0.1, staplewidth = 0.1, color = "black") +
  geom_jitter(aes(color = Point_Color, shape = Point_Shape), width = 0.12, size = 0.8, alpha = 0.9) +
  scale_color_identity() +
  scale_shape_identity() + 
  theme_classic(base_size = 11) +
  scale_fill_manual(values = box_colors) +
  labs(
    title = NULL,
    x = "Genomic Compartment & Sample Source",
    y = "Chao 1 Alpha Diversity (ARGs)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = comparisons_list,
                     method = "wilcox.test",
                     label = "p.signif", 
                     bracket.size = 0.35,
                     tip.length = 0.015,
                     size = 2,
                     step.increase = 0.07) +
  coord_cartesian(ylim = c(0, max_y * 1.7))

# Save high-resolution plot to PDF
pdf(OUTPUT_PLOT, width = 2.5, height = 4)
print(p)
dev.off()

cat("✅ Unified 6-box plot successfully saved to:", OUTPUT_PLOT, "\n")
