# Figure 7: Circular Phylogenetic Tree with Metadata Rings

**Description:**
This document outlines the reproduction of Figure 7, a circular phylogenetic tree (nested pie-chart aesthetic) visualizing the Achtman 7-gene MLSTs, phylogroups, and serotypes of 68 ESBL-*E. coli* isolates (45 study isolates + 23 globally published genomes).

### Required Data Files
Located in the `data_main_figures/Figure_7/` directory:
1. `Fig_7_tree.nwk`: The Newick tree file.
2. `metadata.xlsx`: Contains Supplementary Table 16 with the metadata (Phylogroups, Sample Source, District, ST, and Serotypes).

### Reproducibility Instructions
1. Navigate to the [iTOL web server](https://itol.embl.de/).
2. Upload the `Fig_7_tree.nwk` file to generate the base circular tree.
3. Using the `metadata.xlsx` file, create standard iTOL color-strip annotation files for each data layer.
4. Drag and drop the annotation files onto the tree canvas to create the concentric rings in the following order (from inside to outside):
   * **Tree Nodes (Internal Color):** Phylotypes
   * **Core Labels:** Sample IDs (Isolates from this study are marked with an asterisk `*`)
   * **Ring 1 (Inner):** Sample Type / Source (Sheep, Environment, Human Blood)
   * **Ring 2:** Sampling District (Chitradurga, Tumakuru, or Global locations)
   * **Ring 3:** Sequence Types (STs)
   * **Ring 4 (Outermost):** Identified Serotypes
5. Export the final annotated visualization as a high-resolution PDF/SVG.
