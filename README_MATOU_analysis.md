# MATOU_v1_metaT nitrogen-gene correlation analysis

This repository contains a complete R workflow for the Tara Oceans **MATOU_v1_metaT**
analysis used to examine co-variation among nitrogen transporter and assimilation
homologues.

## What the script does

`MATOU_metatranscriptome_analysis_complete.R` performs the full downstream workflow:

1. Reads one gene-specific MATOU_v1_metaT abundance matrix per gene.
2. Aggregates abundance values at seven taxonomic levels:
   - Eukaryota
   - Stramenopiles
   - Bacillariophyta
   - Bacillariophycidae
   - Coscinodiscophyceae
   - Chaetocerotaceae
   - Thalassiosirales
3. Saves a seven-level abundance summary for every gene.
4. Merges all genes into a single abundance table.
5. Calculates pairwise **Spearman correlations from raw abundance values**.
6. Performs complete-linkage hierarchical clustering using **1 - rho** as the distance.
7. Generates:
   - the seven correlation heatmaps corresponding to Supplementary Fig. S11a;
   - the Eukaryota, Bacillariophyta and Chaetocerotaceae heatmaps used in Fig. 4a;
   - the focused UT1-NR Spearman correlation plots corresponding to Supplementary Fig. S11b.
8. Saves exact Spearman rho, P values, pairwise n values and R session information.

## Important statistical point

Spearman correlations are calculated from the **original, untransformed abundance values**.

For the UT1-NR scatter plots only, points are displayed after:

```text
log10(abundance + 1 × 10^-10)
```

The fitted straight line and confidence band on those plots are a visual guide on the
log-transformed coordinates and are not used to calculate Spearman rho.

## Upstream homologue identification

The script starts from gene-specific MATOU abundance matrices. It does **not** reproduce
the upstream Ocean Gene Atlas / BLASTP search.

For the manuscript analysis, homologues were identified in MATOU_v1_metaT using BLASTP
with an E-value threshold of `1 × 10^-10`, and the resulting abundance matrices were used
as inputs for this R workflow.

## Required R packages

```r
install.packages(c(
  "readxl",
  "openxlsx",
  "dplyr",
  "tidyr",
  "ggplot2",
  "pheatmap",
  "gridExtra"
))
```

## Recommended directory structure

```text
project/
├── MATOU_metatranscriptome_analysis_complete.R
├── data/
│   └── raw/
│       ├── abundance_matrix_UT1.xlsx
│       ├── abundance_matrix_UT2.xlsx
│       ├── abundance_matrix_UT3.xlsx
│       ├── abundance_matrix_UT4.xlsx
│       ├── abundance_matrix_NR.xlsx
│       └── ... one file per gene
└── results/
```

The `results/` directory is created automatically if it does not exist.

## Input format

Each input Excel workbook should contain one abundance matrix with:

- column 1: `Gene_ID`
- column 2: `taxonomy`
- columns 3 onward: MATOU_v1_metaT sample abundance values

The uploaded UT3 example follows this structure.

## Gene naming

The script tries to infer a gene name from either:

- the filename (for example `abundance_matrix_UT3.xlsx`), or
- the first worksheet name (for example `abundance_matrix-UT3.`).

If a filename or worksheet does not map cleanly to the manuscript gene name, add a
manual entry in `gene_name_overrides` near the top of the script.

Example:

```r
gene_name_overrides <- c(
  "abundance_matrix_UreDF.xlsx" = "UreD/F"
)
```

## Main output folders

```text
results/
├── 01_taxonomy_aggregated/
├── 02_correlation_matrices/
└── 03_figures/
```

Notable files include:

- `all_genes_taxonomy_merged_7levels.xlsx`
- `Spearman_rho_<taxon>.csv`
- `Spearman_pvalue_<taxon>.csv`
- `Spearman_n_<taxon>.csv`
- `UT1_NR_Spearman_statistics.csv`
- `Fig4a_representative_Spearman_heatmaps.pdf`
- `Supplementary_Fig_S11a_pairwise_Spearman_heatmaps.pdf`
- `Supplementary_Fig_S11b_UT1_NR_Spearman_scatter.pdf`
- `Supplementary_Fig_S11_complete_layout.pdf`
- `analysis_sessionInfo.txt`

## Running the analysis

Set the R working directory to the project root and run:

```r
source("MATOU_metatranscriptome_analysis_complete.R")
```

or from a terminal with R installed:

```bash
Rscript MATOU_metatranscriptome_analysis_complete.R
```

## Notes on figure appearance

The script preserves the analysis logic and generates publication-ready versions of the
correlation heatmaps and UT1-NR plots. Final panel dimensions, typography and assembly
may still be adjusted in a graphics editor for journal layout without changing the data
analysis.
