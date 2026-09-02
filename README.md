# NR–Urea analysis

Code and data-analysis workflows for the study:

**“Nitrate reductase-dependent control of urea acquisition in diatoms”**

This repository contains custom scripts used for RNA-seq and Tara Oceans
MATOU_v1_metaT metatranscriptomic analyses.

## Contents

### RNA-seq analysis
Scripts for:
- HISAT2 read alignment
- samtools processing
- featureCounts quantification
- DESeq2 differential-expression analysis
- KEGG and GO enrichment analysis

Raw RNA-seq data are available from the NCBI Sequence Read Archive under
BioProject **PRJNA1420504**.

### MATOU_v1_metaT analysis
Scripts for:
- taxonomic aggregation of gene-specific abundance matrices
- pairwise Spearman correlation analysis
- correlation heatmaps used in Fig. 4a and Supplementary Fig. S11
- focused UT1–NR correlation analysis
