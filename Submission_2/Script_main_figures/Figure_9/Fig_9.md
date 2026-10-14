# Figure 9: Gene Neighbourhood and Synteny Analysis

**Description:**
This directory contains interactive HTML visualizations (Figures 9a–9e) showing the genomic context and synteny alignments for *blaCTX-M-15*, Type IV pilus, iron uptake, and toxin virulence genes. These files can be opened and viewed interactively in any standard web browser.

Because these visualizations were generated directly from our bioinformatics pipeline (bakta followed by clinker), the track labels contain embedded metadata. Please use the guide below to interpret the labels for each figure.

### 1. Naming Convention for Figures 9a & 9b
*(10 kb upstream/downstream synteny of blaCTX-M-15 on plasmid and chromosomal contigs)*

The track labels in these HTML files follow a structured format containing all relevant metadata for that specific contig. 
**Example Label:** `RS_S48_blaCTX-M-15_2290020_1706_2581_BA_17`
* **`RS`**: Sample Type (e.g., RS = Rectal Swab)
* **`S48`**: Internal Sequencing Pipeline ID (Seq ID)
* **`blaCTX-M-15`**: Target Gene
* **`2290020`**: Contig ID
* **`1706`**: Start coordinate on the contig
* **`2581`**: End coordinate on the contig
* **`BA`**: Village ID
* **`17`**: Final Sample ID as reported in the manuscript

### 2. Naming Convention for Figures 9c, 9d & 9e
*(Arrangement of Type IV pilus, iron uptake, and toxin virulence genes)*

The track labels in these HTML files combine the target gene cluster and the internal sequencing ID.
**Example Labels:** `cofS73_plasmid.fasta`, `iucsitS34_plasmid.fasta`, `hlyS58_plasmid.fasta`
* **`cof` / `iucsit` / `hly`**: The target gene or gene cluster (e.g., *cof* = Type IV pilus, *iuc/sit* = iron uptake, *hly* = toxin).
* **`S73` / `S34` / `S58`**: Internal Sequencing Pipeline ID (Seq ID).
* **`_plasmid.fasta`**: Denotes the genomic compartment of the assembly.

**Mapping to Manuscript IDs (Crucial for Fig 9c, 9d, 9e):**
The IDs displayed in these specific alignments (e.g., `S73`) are the internal Seq IDs. To map these to the final **Manuscript Sample IDs**, please cross-reference the provided `Metadata.xlsx` file. The metadata table provides a direct translation from `Sample ID[SEQ_FILE]` to `Sample ID [AS PER MANUSCRIPT]`.
