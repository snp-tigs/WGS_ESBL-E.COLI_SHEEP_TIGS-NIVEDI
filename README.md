# **Whole-genome analysis reveals shared extended-spectrum β-lactamase-producing *Escherichia coli* lineages between migratory sheep and surrounding farm environments** <img src="/Graphical_Abstract.png" alt="Graphical Abstract" align="center" width="1200"> 
## **Abstract**
The spillover of pathogenic organisms across livestock-environmental interface is a major One Health concern. To investigate the co-occurrence, antimicrobial resistance profiles, and genomic characteristics of extended-spectrum β-lactamase-producing _Escherichia coli_ (ESBL-_E. coli_) in healthy migratory sheep flocks and associated farm environments, a cross-sectional study was conducted in South Karnataka, India. 9124 animals from 60 farms were included in the study, and 701 samples from sheep and environment were analysed, during 2023-2025. The occurrence of ESBL-_E. coli_ was comparable in sheep and environmental samples, and exhibited multidrug resistance patterns. Whole-genome-sequencing of ESBL-_E. coli_ isolates revealed substantial diversity in antimicrobial resistance genes on chromosomes and plasmids, including efflux pumps and a wide range of β-lactamase genes, including _blaCTX-M, blaTEM, blaDHA, AmpC_, and _blaSHV_ variants. Several plasmid replicon types were detected in both sheep and environmental isolates, suggesting potential plasmid-mediated horizontal gene transfer at the livestock-environment interface. Overall, the detection of multidrug-resistant ESBL-_E. coli_ in migratory sheep and their surrounding environments, together with shared resistomes, plasmids, and lineages, underscores the role of small non-sedentary ruminants as overlooked nodes in AMR dissemination and highlights the need for targeted One Health surveillance strategies. The study indicates the sharing of E. coli lineages across the animals- farm interface.


## **Overview**
This repository contains the datasets, metadata, and computational scripts required to reproduce the genomic analyses presented in our manuscript. The repository provides full transparency for peer review, allowing for the independent reproduction of all phenotypic heatmaps, diversity indices, geospatial maps, phylogenetic trees, and genomic synteny alignments.

## **📂 Repository Structure**

The repository is organized into two primary directories reflecting the peer-review timeline:

### **1. `Submission_1/`**
Contains the original raw datasets, initial scripts for figures, and datasets provided during the first round of manuscript submission. 

### **2. `Submission_2/` (Current Revision - Oct 2026)**
Contains the highly standardized, finalized datasets and refactored computational scripts generated for the revised manuscript. This folder represents the most up-to-date and reproducible version of the project.

Inside `Submission_2/`, the workflow is strictly divided for clarity:
- **`data_main_figures/`**: Contains all raw inputs organized by figure number (Figures 1–9), including CSV presence/absence matrices, GeoJSON boundaries, Newick (`.nwk`) trees, and master metadata files.
- **`script_main_figures/`**: Contains the execution scripts (`.R` and `.py`) and methodological instructions (`.md`) to generate the publication-ready visualizations (PDFs, SVGs, and interactive HTMLs). 

## **💻 Usage Instructions**
To reproduce the updated figures from `Submission_2`:
1. Clone this repository to your local machine.
2. Set your working directory to the root of the `Submission_2/` folder.
3. Execute the R or Python scripts directly from the `script_main_figures/` directory. All scripts are programmed using relative paths and will automatically export high-resolution vector graphics to locally generated `results/` folders.

## **Contact**  

| Name | Role | Affiliation | Contact |
| :--- | :--- | :--- | :--- |
| **Shivranjani C Moharir, Ph.D.** | Senior Scientist | Tata Institute for Genetics and Society (TIGS), Bangalore | 📧 [shivranjani.c@tigs.res.in](mailto:shivranjani.c@tigs.res.in) |
| **Manas Madhukar** | Ph.D. Student | Tata Institute for Genetics and Society (TIGS), Bangalore | 📧 [manas.madhukar@tigs.res.in](mailto:manas.madhukar@tigs.res.in) |
