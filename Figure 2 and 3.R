library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(RColorBrewer)
library(grid)

#laoding data
df <- read_delim("/mnt/406A19506A1943D6/TIGS_NIVEDI/Main Figures_4_Feb_2026/Figure_2_3/wgs_72_Figure_2_3_Final_10022026.tsv", delim = "\t")

#gene groups
amr_genes <- c(	"Amoxicillin-Clavulanate","Ampicillin","Aztreonam","Cefotaxime","Cefoxitin","Ceftriaxone","Ceftazidime","Imipenam","Choloramphinicol","Tetracycline","Nalidixic acid ", "Enrofloxacin", "Co-trimoxazole","Amikacin")
esbl_genes <- c("TEM", "SHV", "AmpC", "CTXM")
virulence_genes <- c("stx1", "stx2")
tet_genes <- c("tetA", "tetB")

#cleaning missing values
gene_cols <- c(amr_genes, esbl_genes, virulence_genes, tet_genes)
df[gene_cols] <- lapply(df[gene_cols], function(x) as.numeric(ifelse(is.na(x) | x == "", 0, x)))

#Sample ID as rownames
df <- column_to_rownames(df, var = "Sample ID")

#convert annotations to factor
df$District <- as.factor(df$District)
df$Village <- as.factor(df$Village)
df$`Sample Type` <- factor(df$`Sample Type`, levels = c( "Rectal swab", "Soil", "Feed", "Water"))  # Explicit order

#Custom color palettes

district_col <-  c(
  "Chitradurga" = "#ff6d3e",
  "Tumkur"      = "#252d80"
)

village_col <- c(
  "Dodderi" = "#c8755b", "Ganjigunte" = "#e18466", "Pavgada ward no 1" = "#636ede",
  "Bharamasagara" = "#fa9372", "Lokkamanahalli" = "#6f7bf7", "Talakere" = "#7d88f7",
  "Maddakkanahalli" = "#8b95f8", "Pavgada ward no 2" = "#9aa2f9", "Nerlegudda" = "#a8affa",
  "Anekere" = "#b7bdfb", "Herur" = "#c5cafb", "Dandikere" = "#d3d7fc",
  "Kilarahalli" = "#434ca9", "Bagganandu Kaval" = "#fa9372", "Doddamalligere" = "#2f39a0",
  "Muddapura" = "#fa9d80", "Shivapura" = "#8288c6", "Thirumalapalya" = "#5862c5",
  "Tumkur ward no-9" = "#2a3390", "Sondekere" = "#fbb39c", "Janapanahalli" = "#2f39a0",
  "Kallankere" = "#8288c6", "Mallapura" = "#fcc9b8", "Hosahalli G" = "#abafd9", "Vanivilasapura" = "#fdded4"
)

stype_col <- c(
  "Feed"         = "#adfda2",
  "Soil"         = "#beb15b",
  "Rectal swab"  = "#5c3863",
  "Water"        = "#4d8eff"
)

#sorting the data frame by Sample Type
df <- df %>%
  arrange(`Sample Type`)

#Row annotations with Sample ID just after Sample Type
metadata_anno <- rowAnnotation(
  District = df$District,
  Village = df$Village,
  SampleType = df$`Sample Type`,
  col = list(
    District = district_col,
    Village = village_col,
    SampleType = stype_col
  ),
  annotation_name_side = "top",
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    District = list(title = "District"),
    Village = list(title = "Village"),
    SampleType = list(title = "Sample Type")
  )
)

sampleid_anno <- rowAnnotation(
  SampleID = anno_text(
    rownames(df),
    gp = gpar(fontsize = 9),
    just = "left"
  ),
  width = max_text_width(rownames(df)) + unit(0.00000001, "mm"),
  show_annotation_name = TRUE
)

#color functions for each heatmap
col_fun_amr <- colorRamp2(c(0, 1), c("white", "#D73027"))
col_fun_esbl <- colorRamp2(c(0, 1), c("white", "#4575B4"))
col_fun_vir <- colorRamp2(c(0, 1), c("white", "#66BD63"))
col_fun_tet <- colorRamp2(c(0, 1), c("white", "#F46D43"))

# heatmaps
hm_amr <- Heatmap(
  as.matrix(df[, amr_genes]),
  name = "Antibiotic Classes",
  col = col_fun_amr,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_title = "AMR genes",
  column_title_side = "top",
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.2),
  column_names_side = "top"
)

hm_esbl <- Heatmap(
  as.matrix(df[, esbl_genes]),
  name = "ESBL genes",
  col = col_fun_esbl,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_title = "ESBL genes",
  column_title_side = "top",
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.2),
  column_names_side = "top"
)

hm_vir <- Heatmap(
  as.matrix(df[, virulence_genes]),
  name = "Virulence genes",
  col = col_fun_vir,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_title = "Virulence genes",
  column_title_side = "top",
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.2),
  column_names_side = "top"
)

hm_tet <- Heatmap(
  as.matrix(df[, tet_genes]),
  name = "Tetracycline genes",
  col = col_fun_tet,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_title = "Tetracycline genes",
  column_title_side = "top",
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.2),
  column_names_side = "top"
)

# complete layout
final_plot <- metadata_anno + 
  sampleid_anno + 
  hm_amr + 
  hm_esbl + 
  hm_vir + 
  hm_tet

draw(
  final_plot,
  heatmap_legend_side = "top",
  annotation_legend_side = "bottom"
)

