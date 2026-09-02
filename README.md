# Chaetoceros muelleri RNA-seq analysis

This directory contains cleaned, reproducible scripts reconstructed from the
original RNA-seq alignment, DESeq2 differential-expression, and KEGG/GO
enrichment scripts used for the study.

## Raw sequencing data

The raw RNA-seq data are deposited in the **NCBI Sequence Read Archive (SRA)**
under BioProject:

**PRJNA1420504**

The scripts do not hard-code individual SRA Run accessions because the uploaded
historical workflow used laboratory/service-provider FASTQ basenames. The file
`sample_metadata_PRJNA1420504.csv` preserves those historical basenames and the
mapping to the biological conditions used in the analysis.

When starting directly from SRA downloads, use the BioProject run table to map
each SRR accession to the corresponding biological sample, then rename or
symlink the paired FASTQ files to the `raw_basename` names in the metadata file.

## Experimental sample structure recovered from the historical script

The alignment script contained:

- **42 biological conditions**
- **3 biological replicates per condition**
- **126 paired-end RNA-seq libraries in total**

The cleaned metadata file contains one row per biological replicate.

Legacy analysis labels are preserved so that the cleaned code remains traceable
to the historical analysis. The metadata also includes the manuscript strain
labels used for the principal strains:

| Legacy label | Manuscript label |
|---|---|
| WT | WT |
| 13509 | 13-NR509 |
| C5 | C5 |
| NNR | NRpNR-1 |
| NUT1 | NRpUT1-12 |
| NUT2 | NRpUT2-5 |

Treatment codes:

| Code | Meaning |
|---|---|
| NO | nitrate |
| UR | urea |
| NH | ammonium |
| NF | nitrogen-free |
| 0H | pre-transfer baseline |

## Files

### `sample_metadata_PRJNA1420504.csv`

Contains:

- clean sample ID;
- historical condition name;
- legacy strain label;
- manuscript strain label;
- nitrogen treatment;
- sampling time;
- biological replicate number;
- historical FASTQ basename;
- BioProject accession.

### `01_align_and_featurecounts.sh`

Reconstructs the upstream read-processing workflow from the uploaded shell
script:

1. optional `fastp` preprocessing;
2. HISAT2 alignment;
3. `samtools sort`;
4. BAM indexing;
5. `samtools flagstat`;
6. optional per-base depth;
7. `featureCounts` quantification.

The historical script used the following featureCounts settings:

```text
-p -C -F GFF -t mRNA -g ID
```

These options are preserved.

### Important note about fastp

The historical shell script contained a `fastp` command with the following
settings:

```text
-q 10 -u 50 -y -g -Y 10 -e 20 -l 100 -b 150 -B 150
```

However, that entire block was **commented out** in the uploaded source file.
Therefore, the cleaned workflow keeps fastp as an optional step and sets:

```bash
RUN_FASTP=false
```

by default. Do not claim that this fastp step was used unless that can be
confirmed from the original preprocessing records.

### `02_deseq2_differential_expression.R`

Performs DESeq2 differential-expression analysis with:

```text
design = ~ condition
```

The two low-expression filters found in the historical script are retained:

```text
rowMeans(count) > 1
count >= 5 in at least 3 samples
```

The script automatically generates two types of contrasts:

1. each treatment/time point versus the corresponding 0-h baseline within the
   same strain;
2. each non-WT strain versus WT under the matched treatment and time point.

This reproduces the contrast logic written out manually in the original script.

## Differential-expression criteria

The historical workflow used two related but distinct thresholds:

- **Overall DEG set:** `P < 0.05` and `|log2FC| >= 1`.
- **Upregulated subset:** `padj < 0.05` and `log2FC >= 1`.
- **Downregulated subset:** `padj < 0.05` and `log2FC <= -1`.

The directional UP/DOWN sets are therefore intentionally more stringent subsets
used when genes were separated by direction, including downstream enrichment
analyses.

The cleaned script writes these results to:

```text
results/deseq2/DEG_ALL_P005_LFC1/
results/deseq2/DEG_UP_PADJ005_LFC1/
results/deseq2/DEG_DOWN_PADJ005_LFC1/
```

### `03_kegg_go_enrichment.R`

Runs custom over-representation analysis with:

```r
clusterProfiler::enricher()
```

and uses `enrichplot::dotplot()` for enrichment visualization.

and preserves the principal settings present in the historical script:

```text
pAdjustMethod = "BH"
maxGSSize = 500
showCategory = 20
```

The script saves the complete enrichment output and also a
BH-adjusted-P `< 0.05` subset.

This step requires gene-to-term mapping tables that are **not contained in
BioProject PRJNA1420504** and therefore need to be supplied separately:

```text
annotations/KEGG_N.tsv
annotations/GO_N.tsv
```

Expected KEGG columns:

```text
ko_id    gene_id    pathway
```

Expected GO columns:

```text
go_id    gene_id    annotation
```

(`pathway` is also accepted instead of `annotation` for the GO term-name
column.)

## Recommended repository layout

```text
RNAseq_analysis_PRJNA1420504/
├── README.md
├── sample_metadata_PRJNA1420504.csv
├── 01_align_and_featurecounts.sh
├── 02_deseq2_differential_expression.R
├── 03_kegg_go_enrichment.R
├── data/
│   └── raw_fastq/
├── reference/
│   ├── Chaetoceros.t2t.fa
│   └── Chaetoceros.gene.as.gff
├── annotations/
│   ├── KEGG_N.tsv
│   └── GO_N.tsv
└── results/
```

Large FASTQ, BAM, and reference files do not need to be committed to GitHub.
The raw reads should instead be referenced by the stable SRA BioProject
accession **PRJNA1420504**.

## Running

### 1. Alignment and counting

```bash
bash 01_align_and_featurecounts.sh
```

Environment variables can override paths, for example:

```bash
RAW_DIR=/path/to/fastq \
REF_FASTA=/path/to/Chaetoceros.t2t.fa \
GFF=/path/to/Chaetoceros.gene.as.gff \
bash 01_align_and_featurecounts.sh
```

### 2. DESeq2

From the project root:

```bash
Rscript 02_deseq2_differential_expression.R
```

### 3. KEGG/GO enrichment

After placing the annotation mappings in `annotations/`:

```bash
Rscript 03_kegg_go_enrichment.R
```

## Reproducibility outputs

The scripts save:

- software versions for HISAT2, samtools and featureCounts;
- the exact sample metadata used;
- DESeq2 contrast manifests;
- filtering summaries;
- normalized counts;
- complete DESeq2 result tables;
- both historical DEG threshold sets;
- full and filtered KEGG/GO enrichment tables;
- R `sessionInfo()` files.

## Files intentionally not used as the public analysis code

The uploaded `.Rhistory` and `.RData` files are session-state/history files.
They are useful for reconstructing the historical workflow, but they are not a
good public reproducibility format. The cleaned scripts above should be
uploaded instead.

## Verification before publication

Before depositing the code on GitHub/Zenodo, verify:

1. the mapping between SRA Run accessions and the `raw_basename` entries;
2. whether the optional fastp step was actually used;
3. the exact HISAT2, samtools, Subread/featureCounts and R package versions;
4. that the reference FASTA/GFF files correspond exactly to those used for the
   published count matrix;
5. that the KEGG/GO mapping tables correspond to the annotation version used in
   the published analysis.
