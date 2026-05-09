library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(tibble)

# data
df <- read.csv("/path/to/input_file/x.csv", check.names = FALSE)

#Extract metadata and matrix
metadata <- df[, 1:4]
amr_matrix <- as.matrix(df[, 5:ncol(df)])
rownames(amr_matrix) <- metadata$`Sample ID`

#custom color palettes
district_col <-  c(
  "Chitradurga" = "#ff6d3e", #af664f",
  "Tumkur"      = "#252d80"
  
)

village_col <- c(
  "Dodderi"        = "#c8755b",
  "Ganjigunte"     = "#e18466",
  "Pavgada ward no 1" = "#636ede", #T
  "Bharamasagara"  = "#fa9372",
  "Lokkamanahalli" = "#6f7bf7", #T
  "Talakere" = "#7d88f7", #T
  "Maddakkanahalli" = "#8b95f8",#T
  "Pavgada ward no 2" = "#9aa2f9",#T
  "Nerlegudda" = "#a8affa", #T
  "Anekere" = "#b7bdfb", #T
  "Herur" = "#c5cafb", #T
  "Dandikere" = "#d3d7fc", #T
  "Kilarahalli" = "#434ca9", #T
  "Bagganandu Kaval" = "#fa9372",
  "Doddamalligere" = "#2f39a0", #T
  "Muddapura" = "#fa9d80",
  "Shivapura" = "#8288c6", #T
  "Thirumalapalya" = "#5862c5",#T
  "Tumkur ward no-9" = "#2a3390", #T
  "Sondekere" = "#fbb39c",
  "Janapanahalli" = "#2f39a0", #T
  "Kallankere" = "#8288c6", #T
  "Mallapura" = "#fcc9b8",
  "Hosahalli G" = "#abafd9", #T
  "Vanivilasapura" = "#fdded4"
)

stype_col <- c(
  "Feed"         = "#adfda2",
  "Soil"         = "#beb15b",
  "Rectal swab"  = "#5c3863",
  "Water"        = "#4d8eff"
)

#Create row annotations
row_annot <- rowAnnotation(
  District = metadata$District,
  Village = metadata$Village,
  SampleType = metadata$`Sample Type`,
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

#manual color bins for heatmap values
col_fun_discrete <- function(x) {
  if (x == 0) return("#D3D3D3")
  else if (x <= 1) return("#5fad8b")
  else if (x <= 5) return("#1b8a5a")        
  else if (x <= 10) return("#fcc763")       
  else if (x <= 15) return("#f68838")
  else if (x <= 20) return("#ee3e32")
  else return("#67000D")                    
}

# Vectorize for matrix input
col_discrete_vectorized <- Vectorize(col_fun_discrete)

#Legend for custom colors
custom_legend <- Legend(
  labels = c("0","0.1-1", "2-5", "6-10", "11-15", "16-20"),
  title = "Normalised Virulence gene copy per genome",
  legend_gp = gpar(fill = c("#D3D3D3","#5fad8b", "#1b8a5a",'#fcc763', "#f68838", "#ee3e32","#67000D"))
)
#heatmap
ht <- Heatmap(
  amr_matrix,
  name = "Normalised RPM",
  col = col_discrete_vectorized,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  row_names_side = "left",
  left_annotation = row_annot,
  column_dend_side = "top",
  row_dend_side = "left",
  rect_gp = gpar(col = "black", lwd = 0.4)
)

#heatmap + custom legend
draw(ht,
     annotation_legend_list = list(custom_legend),
     heatmap_legend_side = "right",
     annotation_legend_side = "right"
)

#Extract Clustered Column Order

#Use column_order() to get the index of columns as they appear in the plot
col_indices <- column_order(ht)

#Extract the actual column names using those indices
clustered_col_names <- colnames(amr_matrix)[col_indices]

# output path
output_txt_path <- "/path/to/output_folder/Clustered_Column_Names.txt"

# optional: Write to a text file (one name per line)
writeLines(clustered_col_names, output_txt_path)

# Optional: Print to console to verify
cat("Success! Clustered column names saved to:", output_txt_path, "\n")
print(head(clustered_col_names))

