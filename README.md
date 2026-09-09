# Project Overview

*Saccharomyces cerevisiae* exhibits notable phenotypic diversity when exposed to environmental stressors such as metal toxicity. This study investigates how divergence in Copy Number Variation (CNV) and downstream expression across five key copper stress regulation genes (*CUP1*, *SOD1*, *HSP12*, *CTR1*, *CRS5*) explains variations in transcriptional plasticity and growth fitness across four diverse *S. cerevisiae* isolates (**Y12**, **YJM1389**, **YPS128**, **YPS163**).
By synthesizing high-throughput RNA-seq data, CNV matrices, and phenotypic fitness profiles, this project demonstrates that strain-specific copper tolerance is mediated specifically through gene-dosage amplification of metallothioneins (*CUP1* and *CRS5*), rather than regulatory divergence in primary transport (*CTR1*) or oxidative stress effectors (*SOD1*, *HSP12*).


## Built With

### Upstream HPC & Alignment Pipeline
- **Bash / Slurm:** Command-line automation and cluster execution
- **SRA Toolkit (`fasterq-dump`):** High-throughput FASTQ extraction
- **Trim Galore / Cutadapt:** Quality trimming (Phred > 20, length > 36bp) and adapter removal
- **STAR Aligner:** Fast RNA-seq spliced alignment to the *S. cerevisiae* S288C reference genome
- **Subread (`featureCounts`):** Genomic feature quantification mapping reads to CDS regions

### Downstream R / Bioconductor Environment
- **R Environment:** Statistical programming, multi-omic integration, and plotting
- **`data.table`:** High-performance list-binding and count matrix assembly (`do.call(cbind, ...)`)
- **`tidyverse` & `ggplot2`:** Reshaping (`pivot_longer`), custom linear regression models (`geom_smooth(method = "lm")`), and multi-panel displays
- **`pheatmap`:** Row-scaled hierarchical clustering with phenotype annotation sidebars
- **`ggfortify`:** Principal Component Analysis (PCA) biplots with strain vectors via `autoplot()`
- **`RColorBrewer`:** Qualitative palette mapping (`Set2`) for strain expression profiling

## Key Findings

- **Primary Fitness Drivers:** Differential expression of metallothioneins (*CUP1* and *CRS5*) account for ~78.87% of the total biological variance (PC1) across strains.
- **Conservation of Effectors:** Expression of *CTR1*, *HSP12*, and *SOD1* exhibited near-zero variance across strains (normalized expression ~1.0), indicating that basal transport and general oxidative mechanisms are highly conserved.
- **Strain Clustering:** Row-scaled hierarchical clustering separated isolates into a high-tolerance (**ACK / Y12** and **AKN / YPS128**) and low-tolerance (**SACE_YCK / YJM1389** and **AVI / YPS163**) phenotypic groups.
- **Secondary Regulatory Pathways:** Despite near-identical metallothionein expression levels, strain **ACK** achieved significantly higher growth in 10 mM $\text{CuSO}_4$ (1.1) compared to **AKN** (0.2), pointing to auxiliary metabolic or regulatory factors contributing to peak fitness under severe toxicity.

## Input Datasets & Acquisition

### 1. RNA-seq Datasets (NCBI SRA)
Raw single-end RNA-sequencing accessions retrieved via SRA Toolkit:
- **Y12 (ACK):** `SRR14156703` - `SRR14156711`
- **YJM1389 (SACE_YCK):** `SRR14156757` - `SRR14156765`
- **YPS128 (AKN):** `SRR14156775` - `SRR14156783`
- **YPS163 (AVI):** `SRR14156784` - `SRR14156792`

### 2. Pangenomic & Phenotypic Matrices
Sourced from the Yeast 1000 Genomes Project:
- **CNV Data:** `genesMatrix_CopyNumber.tab` (ORFs: *CUP1* [`YGR078C-A`], *CRS5* [`YOR031W`], *CTR1* [`YOL112W`], *SOD1* [`YJR104C`], *HSP12* [`YEL024W`])
- **Phenotype Matrix:** `pheno_35Conditions_NormalizedByYPD.txt` (Condition: 10 mM $\text{CuSO}_4$ growth ratio relative to YPD at 30°C)

## Bioinformatic Pipeline

### Phase 1: High-Performance Computing Pipeline (UNIX / HPC Cluster)
Raw FASTQ (SRA) ──► Trim Galore (Phred > 20, >36bp) ──► STAR Alignment (S288C) ──► featureCounts (CDS)

### Phase 2: Statistical Modeling & Visualization (R / Bioconductor)
Step 1: File Path Ingestion & Counts Binding (`list.files()` + `data.table`)

──► Step 2: Subsetting target strains & building unified `master_df`

──► Step 3: Linear regression modeling (`geom_smooth()`) & expression profiling

──► Step 4: Zero-variance filtering & row-scaled `pheatmap` clustering

──► Step 5: Scaled `prcomp()` ordination & `ggfortify` PCA biplots

──► Step 6: *CUP1* vs *CRS5* differentiator analysis ──► Export to `outputs/*.png`
