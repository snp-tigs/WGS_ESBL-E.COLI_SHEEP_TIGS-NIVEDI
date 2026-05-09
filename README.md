# WGS_ESBL-E.COLI_SHEEP_TIGS-NIVEDI
Whole-genome analysis reveals exchange of ESBL-E. coli genomes between migratory sheep and surrounding farm environments

This repository contains the analytical pipeline and visualization scripts used to investigate the occurrence, genomic characteristics, and antimicrobial resistance (AMR) profiles of ESBL-producing E. coli in migratory sheep and their environments in South Karnataka, India.

The code supports a cross-sectional OneHealth study (2023–2025) involving WGS (Whole Genome Sequencing) of isolates to identify shared resistomes and plasmid-mediated gene transfer.

Core Analytical Modules
1. Data Normalization & Quantitation
Normalisation_ARGs_copy_no_per_genome.py: Calculates the abundance of Antimicrobial Resistance Genes (ARGs) by normalizing copy numbers per bacterial genome, allowing for quantitative comparison across diverse sheep and environmental samples.

Normalisation_ARGs_Vir_genes_copy_no_per_genome.py: Extends the normalization logic to include Virulence Factors (VFs), facilitating the study of the "pathobiome" and the public health risks associated with stx-carrying isolates.

2. Diversity & Ecological Analysis
Alpha_diversity_plot.py: Computes and visualizes the within-sample diversity of the resistomes and lineages, helping to compare the microbial complexity between sedentary and migratory livestock environments.

3. Statistical Visualization (Results Reproduction)
Figure 2 and 3.R: R-based scripts for generating core comparative figures, likely focusing on the distribution of ESBL variants (CTX-M, TEM, etc.) and phenotypic resistance patterns across the 60 sampled farms.

Figure 4 to 7.R: Advanced visualization scripts for comparative genomics, plasmid replicon distribution, and phylogenetic lineages, highlighting the horizontal gene transfer at the livestock-environment interface.

4. Geospatial Documentation
Sample_sites_map_graphical_abstract.py: A Python workflow (using geopandas and matplotlib) to generate high-resolution maps of the South Karnataka study area. This script produces the spatial context for the 701 samples collected, specifically designed for the paper's graphical abstract.
