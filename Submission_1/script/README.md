# Whole-genome analysis reveals exchange of ESBL-E. coli genomes between migratory sheep and surrounding farm environments


This repository contains the analytical pipeline and visualization scripts used to investigate the occurrence, genomic characteristics, and antimicrobial resistance (AMR) profiles of ESBL-producing E. coli in migratory sheep and their environments in South Karnataka, India.

The code supports a cross-sectional OneHealth study (2023–2025) involving WGS (Whole Genome Sequencing) of isolates to identify shared resistomes and plasmid-mediated gene transfer.

Core Analytical Modules
1. Data Normalization & Quantitation

Normalisation_ARGs_copy_no_per_genome.py: Calculates the abundance of Antimicrobial Resistance Genes (ARGs), Drug classes, Resistance mechanisms by normalizing copy numbers per bacterial genome

Normalisation_Vir_genes_copy_no_per_genome.py: Calculates the abundance of Virulence genes by normalizing copy numbers per bacterial genome

2. Diversity & Ecological Analysis
   
Alpha_diversity_plot.py: visualizes Alpha diversity of ARGs, Drug classes and Virulence genes

4. Statistical Visualization (Results Reproduction)

Figure 2 and 3.R: R-based scripts for generating figure 2 and 3.

Figure 4 to 7.R: R-based scripts for genetation figures 4, 5, 6 and 7

4. Geospatial Documentation

Sample_sites_Figure_1.py: A Python  script to generate high-resolution maps (Figure 1) of the South Karnataka study area which includes sampling coordinates in Tumakuru and Chitradurga. 
