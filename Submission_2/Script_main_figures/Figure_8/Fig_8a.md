# Figure 8a: Maximum Likelihood Phylogenetic Tree

**Description:** 
This document outlines the steps to reproduce the core genome Maximum Likelihood phylogenetic tree (Figure 8a) of the 72 ESBL *E. coli* isolates, rooted against the *Escherichia coli* str. K-12 substr. MG1655 reference genome. 

The visualization was generated using the Interactive Tree Of Life (iTOL) web server.

### Required Data Files
All required files are located in the `data_main_figures/Figure_8/` directory:
1. `Fig_8a_tree_file.nwk` (The core genome alignment tree in Newick format)
2. `Fig_8a_itol_district_colors.txt` (iTOL annotation dataset for District colors)
3. `Fig_8a_itol_village_colors.txt` (iTOL annotation dataset for Village colors)
4. `Fig_8a_itol_sampletype_colors.txt` (iTOL annotation dataset for Sample Type colors)

### Step-by-Step Reproducibility Instructions

**Step 1: Upload the Tree**
1. Navigate to the [iTOL website](https://itol.embl.de/).
2. Click on **Upload** in the top navigation bar.
3. Drag and drop `Fig_8a_tree_file.nwk` into the upload box.

**Step 2: Apply Annotations**
Once the tree is rendered on your screen, you can automatically apply the metadata color strips by dragging and dropping the annotation files directly onto the tree canvas one by one:
1. Drag and drop `Fig_8a_itol_sampletype_colors.txt`.
2. Drag and drop `Fig_8a_itol_village_colors.txt`.
3. Drag and drop `Fig_8a_itol_district_colors.txt`.

**Step 3: Formatting (Matching the Publication)**
* **Tree Layout:** Ensure the tree is set to "Rectangular" based on the final published aesthetic.
* **Scale:** Enable the Tree Scale to represent the divergence of the isolates from the reference genome.
* **Export:** Click the **Export** tab in the control panel to save the resulting annotated tree as a high-resolution SVG or PDF.
