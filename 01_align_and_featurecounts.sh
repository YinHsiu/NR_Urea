#!/usr/bin/env bash
# ==============================================================================
# 01_align_and_featurecounts.sh
#
# Cleaned RNA-seq alignment and read-counting workflow.
#
# Raw RNA-seq data:
#   NCBI BioProject PRJNA1420504
#
# Historical workflow reconstructed from the authors' original analysis script.
# It uses three biological replicates per condition, HISAT2 alignment, samtools
# BAM processing/QC, and featureCounts quantification.
#
# IMPORTANT:
#   The historical script contained a fastp block, but it was commented out.
#   Therefore fastp is preserved here as OPTIONAL and is disabled by default.
#   Enable it only if these parameters reflect the preprocessing actually used.
#
# Expected project structure:
#
#   project/
#   ├── 01_align_and_featurecounts.sh
#   ├── sample_metadata_PRJNA1420504.csv
#   ├── data/
#   │   └── raw_fastq/
#   │       ├── <raw_basename>_1.fq.gz
#   │       └── <raw_basename>_2.fq.gz
#   ├── reference/
#   │   ├── Chaetoceros.t2t.fa
#   │   └── Chaetoceros.gene.as.gff
#   └── results/
#
# The metadata table preserves the historical raw-file basenames and maps them
# to clean sample IDs.
# ==============================================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

METADATA="${METADATA:-${PROJECT_DIR}/sample_metadata_PRJNA1420504.csv}"
RAW_DIR="${RAW_DIR:-${PROJECT_DIR}/data/raw_fastq}"
REF_FASTA="${REF_FASTA:-${PROJECT_DIR}/reference/Chaetoceros.t2t.fa}"
GFF="${GFF:-${PROJECT_DIR}/reference/Chaetoceros.gene.as.gff}"

RESULT_DIR="${RESULT_DIR:-${PROJECT_DIR}/results/alignment}"
BAM_DIR="${BAM_DIR:-${RESULT_DIR}/bam}"
QC_DIR="${QC_DIR:-${RESULT_DIR}/qc}"
COUNT_DIR="${COUNT_DIR:-${PROJECT_DIR}/results/counts}"
INDEX_DIR="${INDEX_DIR:-${PROJECT_DIR}/reference/hisat2_index}"
INDEX_PREFIX="${INDEX_PREFIX:-${INDEX_DIR}/Chaetoceros_t2t}"

THREADS_ALIGN="${THREADS_ALIGN:-20}"
THREADS_SORT="${THREADS_SORT:-12}"
THREADS_COUNT="${THREADS_COUNT:-10}"
SORT_MEMORY="${SORT_MEMORY:-2G}"

# Optional steps.
RUN_FASTP="${RUN_FASTP:-false}"
RUN_DEPTH="${RUN_DEPTH:-true}"

mkdir -p \
  "${BAM_DIR}" \
  "${QC_DIR}" \
  "${COUNT_DIR}" \
  "${INDEX_DIR}"

for cmd in hisat2 hisat2-build samtools featureCounts; do
  command -v "${cmd}" >/dev/null 2>&1 || {
    echo "ERROR: required command not found: ${cmd}" >&2
    exit 1
  }
done

if [[ "${RUN_FASTP}" == "true" ]]; then
  command -v fastp >/dev/null 2>&1 || {
    echo "ERROR: RUN_FASTP=true but fastp was not found." >&2
    exit 1
  }
fi

[[ -f "${METADATA}" ]] || { echo "ERROR: metadata not found: ${METADATA}" >&2; exit 1; }
[[ -f "${REF_FASTA}" ]] || { echo "ERROR: reference FASTA not found: ${REF_FASTA}" >&2; exit 1; }
[[ -f "${GFF}" ]] || { echo "ERROR: GFF not found: ${GFF}" >&2; exit 1; }

# ------------------------------------------------------------------------------
# Record software versions
# ------------------------------------------------------------------------------

{
  echo "BioProject: PRJNA1420504"
  echo
  echo "HISAT2:"
  hisat2 --version | head -n 1
  echo
  echo "samtools:"
  samtools --version | head -n 2
  echo
  echo "featureCounts:"
  featureCounts -v 2>&1 | head -n 2
  if [[ "${RUN_FASTP}" == "true" ]]; then
    echo
    echo "fastp:"
    fastp --version 2>&1 | head -n 1
  fi
} > "${RESULT_DIR}/software_versions.txt"

# ------------------------------------------------------------------------------
# Build HISAT2 index if it is not already present
# ------------------------------------------------------------------------------

if [[ ! -f "${INDEX_PREFIX}.1.ht2" && ! -f "${INDEX_PREFIX}.1.ht2l" ]]; then
  echo "[INFO] Building HISAT2 index..."
  hisat2-build "${REF_FASTA}" "${INDEX_PREFIX}"
else
  echo "[INFO] Existing HISAT2 index found."
fi

# ------------------------------------------------------------------------------
# Align all samples listed in metadata
# ------------------------------------------------------------------------------

declare -a BAM_FILES=()

# read CSV columns:
# sample_id,condition,strain_legacy,manuscript_strain,treatment_code,treatment,
# time_h,biological_replicate,raw_basename,bioproject
while IFS=',' read -r \
  sample_id condition strain_legacy manuscript_strain treatment_code treatment \
  time_h biological_replicate raw_basename bioproject
do
  [[ "${sample_id}" == "sample_id" ]] && continue
  [[ -z "${sample_id}" ]] && continue

  echo "[INFO] Processing ${sample_id} (${condition}, replicate ${biological_replicate})"

  raw_r1="${RAW_DIR}/${raw_basename}_1.fq.gz"
  raw_r2="${RAW_DIR}/${raw_basename}_2.fq.gz"

  [[ -f "${raw_r1}" ]] || { echo "ERROR: missing FASTQ: ${raw_r1}" >&2; exit 1; }
  [[ -f "${raw_r2}" ]] || { echo "ERROR: missing FASTQ: ${raw_r2}" >&2; exit 1; }

  if [[ "${RUN_FASTP}" == "true" ]]; then
    clean_r1="${QC_DIR}/${sample_id}.clean_1.fq.gz"
    clean_r2="${QC_DIR}/${sample_id}.clean_2.fq.gz"

    # Parameters preserved from the commented fastp block in the original script.
    fastp \
      -q 10 \
      -u 50 \
      -y \
      -g \
      -Y 10 \
      -e 20 \
      -l 100 \
      -b 150 \
      -B 150 \
      -i "${raw_r1}" \
      -I "${raw_r2}" \
      -o "${clean_r1}" \
      -O "${clean_r2}" \
      --html "${QC_DIR}/${sample_id}.fastp.html" \
      --json "${QC_DIR}/${sample_id}.fastp.json"

    input_r1="${clean_r1}"
    input_r2="${clean_r2}"
  else
    input_r1="${raw_r1}"
    input_r2="${raw_r2}"
  fi

  sam_file="${BAM_DIR}/${sample_id}.sam"
  bam_file="${BAM_DIR}/${sample_id}.sorted.bam"

  hisat2 \
    -p "${THREADS_ALIGN}" \
    -x "${INDEX_PREFIX}" \
    -1 "${input_r1}" \
    -2 "${input_r2}" \
    -S "${sam_file}" \
    2> "${QC_DIR}/${sample_id}.hisat2.log"

  samtools sort \
    -@ "${THREADS_SORT}" \
    -m "${SORT_MEMORY}" \
    -o "${bam_file}" \
    "${sam_file}"

  rm -f "${sam_file}"

  samtools index "${bam_file}"

  samtools flagstat \
    "${bam_file}" \
    > "${QC_DIR}/${sample_id}.flagstat.txt"

  if [[ "${RUN_DEPTH}" == "true" ]]; then
    samtools depth -aa \
      "${bam_file}" \
      > "${QC_DIR}/${sample_id}.depth.txt"
  fi

  BAM_FILES+=("${bam_file}")

done < "${METADATA}"

if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
  echo "ERROR: no samples were processed." >&2
  exit 1
fi

# ------------------------------------------------------------------------------
# Quantify reads with featureCounts
# ------------------------------------------------------------------------------

COUNT_OUT="${COUNT_DIR}/Chaetoceros_featureCounts.txt"

echo "[INFO] Running featureCounts on ${#BAM_FILES[@]} BAM files..."

# Options preserved from the historical script:
#   -p          paired-end mode
#   -C          do not count chimeric fragments
#   -F GFF      GFF annotation format
#   -t mRNA     count mRNA features
#   -g ID       use the ID attribute as feature identifier
#
# Note: featureCounts behavior for paired-end counting can depend on Subread
# version. The exact historical command is retained here for provenance.
featureCounts \
  -T "${THREADS_COUNT}" \
  -p \
  -C \
  -F GFF \
  -a "${GFF}" \
  -o "${COUNT_OUT}" \
  -t mRNA \
  -g ID \
  "${BAM_FILES[@]}"

echo "[INFO] Alignment and counting completed."
echo "[INFO] Count matrix: ${COUNT_OUT}"
