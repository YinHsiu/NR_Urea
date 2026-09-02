# ==============================================================================
# 03_kegg_go_enrichment.R
#
# Cleaned over-representation analysis for DESeq2 DEG sets using
# clusterProfiler::enricher().
#
# Raw RNA-seq data underlying the DEG analysis:
#   NCBI BioProject PRJNA1420504
#
# Direction-specific DEG sets entering enrichment are defined as:
#   Upregulated:   padj < 0.05 and log2FoldChange >= 1
#   Downregulated: padj < 0.05 and log2FoldChange <= -1
#
# Required annotation mapping tables:
#
#   annotations/KEGG_N.tsv
#       required columns: ko_id, gene_id, pathway
#
#   annotations/GO_N.tsv
#       required columns: go_id, gene_id, and either annotation or pathway
#
# The historical scripts used:
#   pAdjustMethod = "BH"
#   maxGSSize = 500
#   showCategory = 20
#
# To preserve all enrichment results while making the significance filtering
# transparent, this script:
#   - runs enricher() with pvalueCutoff = 1 and qvalueCutoff = 1;
#   - saves the complete enrichment table;
#   - additionally saves a table filtered at BH-adjusted P < 0.05;
#   - plots up to 20 significant categories where available.
# ==============================================================================

required_packages <- c(
  "clusterProfiler",
  "enrichplot",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})


# ==============================================================================
# 0. CONFIGURATION
# ==============================================================================

# Direction-specific DEG sets used for enrichment in the historical workflow:
#   Upregulated:   padj < 0.05 and log2FoldChange >= 1
#   Downregulated: padj < 0.05 and log2FoldChange <= -1
UP_DEG_DIR <- file.path(
  "results",
  "deseq2",
  "DEG_UP_PADJ005_LFC1"
)

DOWN_DEG_DIR <- file.path(
  "results",
  "deseq2",
  "DEG_DOWN_PADJ005_LFC1"
)

KEGG_FILE <- file.path(
  "annotations",
  "KEGG_N.tsv"
)

GO_FILE <- file.path(
  "annotations",
  "GO_N.tsv"
)

OUTPUT_DIR <- file.path(
  "results",
  "enrichment"
)

PADJ_DISPLAY_CUTOFF <- 0.05
SHOW_CATEGORY <- 20
MAX_GS_SIZE <- 500

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "KEGG"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "GO"), recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. INPUT HELPERS
# ==============================================================================

read_annotation_table <- function(file) {

  if (!file.exists(file)) {
    stop("Annotation file not found: ", file)
  }

  ext <- tolower(
    tools::file_ext(file)
  )

  if (ext %in% c("tsv", "txt")) {

    read.delim(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else if (ext == "csv") {

    read.csv(
      file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else {
    stop(
      "Unsupported annotation format: ",
      ext,
      ". Use TSV/TXT or CSV."
    )
  }
}


validate_columns <- function(
    x,
    required,
    label
) {

  missing <- setdiff(
    required,
    colnames(x)
  )

  if (length(missing) > 0) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}


kegg <- read_annotation_table(
  KEGG_FILE
)

go <- read_annotation_table(
  GO_FILE
)

validate_columns(
  kegg,
  c(
    "ko_id",
    "gene_id",
    "pathway"
  ),
  "KEGG mapping"
)

validate_columns(
  go,
  c(
    "go_id",
    "gene_id"
  ),
  "GO mapping"
)

go_name_column <- if ("annotation" %in% colnames(go)) {
  "annotation"
} else if ("pathway" %in% colnames(go)) {
  "pathway"
} else {
  stop(
    "GO mapping must contain either an 'annotation' or 'pathway' column."
  )
}

KEGG_TERM2GENE <- unique(
  kegg[
    ,
    c(
      "ko_id",
      "gene_id"
    )
  ]
)

KEGG_TERM2NAME <- unique(
  kegg[
    ,
    c(
      "ko_id",
      "pathway"
    )
  ]
)

GO_TERM2GENE <- unique(
  go[
    ,
    c(
      "go_id",
      "gene_id"
    )
  ]
)

GO_TERM2NAME <- unique(
  go[
    ,
    c(
      "go_id",
      go_name_column
    )
  ]
)

colnames(
  GO_TERM2NAME
) <- c(
  "go_id",
  "term_name"
)


# ==============================================================================
# 2. ENRICHMENT AND EXPORT FUNCTIONS
# ==============================================================================

run_enrichment <- function(
    genes,
    term2gene,
    term2name
) {

  genes <- unique(
    genes[
      !is.na(genes) &
        genes != ""
    ]
  )

  if (length(genes) == 0) {
    return(NULL)
  }

  clusterProfiler::enricher(
    gene = genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    maxGSSize = MAX_GS_SIZE
  )
}


export_enrichment <- function(
    enrich_obj,
    output_stem,
    title
) {

  if (is.null(enrich_obj)) {
    writeLines(
      "No input genes were available for enrichment.",
      paste0(
        output_stem,
        "_NO_INPUT_GENES.txt"
      )
    )
    return(
      data.frame(
        n_terms = 0,
        n_significant = 0
      )
    )
  }

  full_result <- as.data.frame(
    enrich_obj@result
  )

  write.csv(
    full_result,
    paste0(
      output_stem,
      "_all.csv"
    ),
    row.names = FALSE
  )

  if (nrow(full_result) == 0) {
    writeLines(
      "No enriched terms were returned.",
      paste0(
        output_stem,
        "_NO_ENRICHED_TERMS.txt"
      )
    )
    return(
      data.frame(
        n_terms = 0,
        n_significant = 0
      )
    )
  }

  significant <- full_result[
    !is.na(full_result$p.adjust) &
      full_result$p.adjust < PADJ_DISPLAY_CUTOFF,
    ,
    drop = FALSE
  ]

  write.csv(
    significant,
    paste0(
      output_stem,
      "_padj_lt_0.05.csv"
    ),
    row.names = FALSE
  )

  if (nrow(significant) > 0) {

    plot_obj <- enrich_obj
    plot_obj@result <- significant

    p <- enrichplot::dotplot(
      plot_obj,
      showCategory = min(
        SHOW_CATEGORY,
        nrow(significant)
      ),
      title = title
    )

    ggsave(
      paste0(
        output_stem,
        "_dotplot.pdf"
      ),
      p,
      width = 4,
      height = 6,
      units = "in"
    )

    ggsave(
      paste0(
        output_stem,
        "_dotplot.png"
      ),
      p,
      width = 4,
      height = 6,
      units = "in",
      dpi = 600
    )
  }

  data.frame(
    n_terms = nrow(full_result),
    n_significant = nrow(significant)
  )
}


# ==============================================================================
# 3. LOOP THROUGH ALL UP/DOWN DEG SETS
# ==============================================================================

up_deg_files <- list.files(
  UP_DEG_DIR,
  pattern = "_UP\\.csv$",
  full.names = TRUE
)

down_deg_files <- list.files(
  DOWN_DEG_DIR,
  pattern = "_DOWN\\.csv$",
  full.names = TRUE
)

deg_files <- c(
  up_deg_files,
  down_deg_files
)

if (length(deg_files) == 0) {
  stop(
    "No directional DEG files found in:\n  ",
    UP_DEG_DIR,
    "\n  ",
    DOWN_DEG_DIR
  )
}

summary_list <- list()

for (deg_file in deg_files) {

  deg <- read.csv(
    deg_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"gene_id" %in% colnames(deg)) {
    stop(
      "DEG file has no gene_id column: ",
      deg_file
    )
  }

  genes <- deg$gene_id

  set_name <- tools::file_path_sans_ext(
    basename(deg_file)
  )

  message(
    "Enrichment: ",
    set_name
  )

  kegg_obj <- run_enrichment(
    genes = genes,
    term2gene = KEGG_TERM2GENE,
    term2name = KEGG_TERM2NAME
  )

  go_obj <- run_enrichment(
    genes = genes,
    term2gene = GO_TERM2GENE,
    term2name = GO_TERM2NAME
  )

  kegg_stats <- export_enrichment(
    kegg_obj,
    file.path(
      OUTPUT_DIR,
      "KEGG",
      paste0(
        set_name,
        "_KEGG"
      )
    ),
    paste0(
      "KEGG: ",
      set_name
    )
  )

  go_stats <- export_enrichment(
    go_obj,
    file.path(
      OUTPUT_DIR,
      "GO",
      paste0(
        set_name,
        "_GO"
      )
    ),
    paste0(
      "GO: ",
      set_name
    )
  )

  summary_list[[length(summary_list) + 1]] <- data.frame(
    DEG_set = set_name,
    n_genes = length(unique(genes)),
    KEGG_terms = kegg_stats$n_terms,
    KEGG_padj_lt_0.05 = kegg_stats$n_significant,
    GO_terms = go_stats$n_terms,
    GO_padj_lt_0.05 = go_stats$n_significant,
    stringsAsFactors = FALSE
  )
}

enrichment_summary <- do.call(
  rbind,
  summary_list
)

write.csv(
  enrichment_summary,
  file.path(
    OUTPUT_DIR,
    "enrichment_summary.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "R_sessionInfo.txt"
  )
)

message(
  "KEGG/GO enrichment completed."
)
