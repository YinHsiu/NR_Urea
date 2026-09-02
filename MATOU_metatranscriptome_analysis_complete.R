# ==============================================================================
# MATOU_metatranscriptome_analysis_complete.R
#
# Complete analysis workflow for nitrogen transporter/assimilation homologues
# in the Tara Oceans MATOU_v1_metaT dataset.
#
# This script reproduces the analysis workflow used for:
#   1) aggregation of each gene-specific MATOU abundance matrix at seven
#      taxonomic levels;
#   2) merging the aggregated abundance matrices across genes;
#   3) pairwise Spearman correlation analysis among all genes;
#   4) correlation heatmaps for the seven taxonomic levels
#      (Supplementary Fig. S11a);
#   5) the three representative heatmaps used in Fig. 4a
#      (Eukaryota, Bacillariophyta, Chaetocerotaceae);
#   6) the focused UT1-NR Spearman correlation scatter plots
#      (Supplementary Fig. S11b).
#
# IMPORTANT STATISTICAL NOTE
# --------------------------
# Spearman correlations are calculated from the ORIGINAL (untransformed)
# MATOU_v1_metaT abundance values.
#
# The transformation
#       log10(abundance + 1e-10)
# is used ONLY for visualization of the UT1-NR scatter plots.
#
# The red line in the scatter plots is a linear fit to the log10-transformed
# values and serves only as a visual guide; it is not the Spearman model.
#
# Upstream homologue identification (e.g. BLASTP in Ocean Gene Atlas with
# E-value < 1e-10) is not performed by this script. The gene-specific abundance
# matrices exported from that analysis are the input files.
#
# Required R packages:
#   readxl, openxlsx, dplyr, tidyr, ggplot2, pheatmap, gridExtra
#
# Suggested repository structure:
#
#   project/
#   ├── MATOU_metatranscriptome_analysis_complete.R
#   ├── data/
#   │   └── raw/
#   │       ├── abundance_matrix_UT1.xlsx
#   │       ├── abundance_matrix_UT2.xlsx
#   │       ├── abundance_matrix_UT3.xlsx
#   │       ├── abundance_matrix_NR.xlsx
#   │       └── ... one file per gene
#   └── results/
#
# Each raw Excel file is expected to contain:
#   column 1: Gene_ID
#   column 2: taxonomy
#   columns 3+: MATOU_v1_metaT sample abundance values
#
# Gene names are inferred from the file name or first worksheet name.
# If necessary, use gene_name_overrides in the CONFIGURATION section below.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURATION
# ==============================================================================

raw_data_dir <- file.path("data", "raw")
output_dir   <- "results"

# All .xlsx files in data/raw are treated as individual gene abundance matrices.
# Keep only raw gene-specific abundance matrices in that directory.
input_file_pattern <- "\\.xlsx$"

# Seven taxonomic levels used in the manuscript.
taxonomy_levels <- data.frame(
  Taxonomy_group = c(
    "Eukaryota",
    "Stramenopiles",
    "Bacillariophyta",
    "Bacillariophycidae",
    "Coscinodiscophyceae",
    "Chaetocerotaceae",
    "Thalassiosirales"
  ),
  Keyword = c(
    "Eukaryota",
    "Stramenopiles",
    "Bacillariophyta",
    "Bacillariophycidae",
    "Coscinodiscophyceae",
    "Chaetocerotaceae",
    "Thalassiosirales"
  ),
  stringsAsFactors = FALSE
)

# Three taxonomic levels displayed in Fig. 4a.
fig4a_taxa <- c(
  "Eukaryota",
  "Bacillariophyta",
  "Chaetocerotaceae"
)

# Three taxonomic levels displayed in Supplementary Fig. S11b.
scatter_taxa <- c(
  "Eukaryota",
  "Bacillariophyta",
  "Chaetocerotaceae"
)

# Pseudocount used only for the scatter-plot visualization.
pseudocount <- 1e-10

# Optional manual correction of gene names.
# The names on the LEFT can be either file basenames or inferred names;
# the values on the RIGHT are the final names to use in the analysis.
#
# Example:
# gene_name_overrides <- c(
#   "abundance_matrix_UreDF.xlsx" = "UreD/F",
#   "UreDF" = "UreD/F"
# )
gene_name_overrides <- c()


# ==============================================================================
# 1. PACKAGE CHECK
# ==============================================================================

required_packages <- c(
  "readxl",
  "openxlsx",
  "dplyr",
  "tidyr",
  "ggplot2",
  "pheatmap",
  "gridExtra"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the script, for example:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(gridExtra)
})


# ==============================================================================
# 2. DIRECTORY SETUP
# ==============================================================================

if (!dir.exists(raw_data_dir)) {
  stop(
    "Input directory not found: ", raw_data_dir,
    "\nCreate data/raw and place one gene-specific MATOU abundance .xlsx file ",
    "per gene in that directory."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "01_taxonomy_aggregated"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "02_correlation_matrices"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "03_figures"), recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

safe_filename <- function(x) {
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}


infer_gene_name <- function(file) {

  file_base <- basename(file)
  file_stem <- tools::file_path_sans_ext(file_base)

  # First try the file name.
  candidate <- sub(
    "^abundance[_-]?matrix[_-]?",
    "",
    file_stem,
    ignore.case = TRUE
  )

  candidate <- sub(
    "_taxonomy_merged.*$",
    "",
    candidate,
    ignore.case = TRUE
  )

  # If the file name is not informative, try the first worksheet name.
  # This also handles sheets such as "abundance_matrix-UT3.".
  if (
    identical(candidate, file_stem) ||
    grepl("^[0-9a-fA-F-]{20,}$", candidate)
  ) {
    sheet_name <- readxl::excel_sheets(file)[1]

    candidate <- sub(
      "^abundance[_-]?matrix[_-]?",
      "",
      sheet_name,
      ignore.case = TRUE
    )
  }

  candidate <- sub("\\.+$", "", candidate)
  candidate <- trimws(candidate)

  # A few common filename-safe forms can be converted back to manuscript names.
  standard_name_map <- c(
    "UreD_F" = "UreD/F",
    "UreD-F" = "UreD/F",
    "Fd_NiR" = "Fd-NiR",
    "NADPH_NiR" = "NADPH-NiR"
  )

  if (candidate %in% names(standard_name_map)) {
    candidate <- unname(standard_name_map[candidate])
  }

  # Optional user-defined overrides.
  if (file_base %in% names(gene_name_overrides)) {
    candidate <- unname(gene_name_overrides[file_base])
  } else if (candidate %in% names(gene_name_overrides)) {
    candidate <- unname(gene_name_overrides[candidate])
  }

  if (is.na(candidate) || candidate == "") {
    stop("Could not infer a gene name from file: ", file)
  }

  candidate
}


taxonomy_match <- function(taxonomy_vector, keyword) {
  # Taxonomy strings are semicolon-delimited. Fixed substring matching mirrors
  # the original analysis while remaining robust to variable lineage depth.
  !is.na(taxonomy_vector) &
    grepl(keyword, taxonomy_vector, fixed = TRUE)
}


aggregate_one_gene <- function(file, gene_name) {

  message("Reading: ", basename(file), "  -> gene: ", gene_name)

  df <- readxl::read_excel(
    file,
    sheet = 1,
    .name_repair = "minimal"
  )

  if (ncol(df) < 3) {
    stop(
      "Input file must contain Gene_ID, taxonomy and at least one abundance ",
      "column: ", file
    )
  }

  colnames(df)[1:2] <- c("Gene_ID", "taxonomy")

  taxonomy <- as.character(df$taxonomy)
  sample_cols <- colnames(df)[3:ncol(df)]

  # Convert abundance columns to numeric.
  abund_num <- as.data.frame(
    lapply(df[, sample_cols, drop = FALSE], function(x) {
      suppressWarnings(as.numeric(as.character(x)))
    }),
    check.names = FALSE
  )

  colnames(abund_num) <- sample_cols

  aggregated_rows <- lapply(
    seq_len(nrow(taxonomy_levels)),
    function(i) {

      group_name <- taxonomy_levels$Taxonomy_group[i]
      keyword    <- taxonomy_levels$Keyword[i]

      idx <- taxonomy_match(taxonomy, keyword)

      # If no homologues match a taxonomic level, use zeros for that gene/level.
      if (sum(idx, na.rm = TRUE) == 0) {
        summed <- rep(0, length(sample_cols))
        names(summed) <- sample_cols
      } else {
        summed <- colSums(
          abund_num[idx, , drop = FALSE],
          na.rm = TRUE
        )
      }

      data.frame(
        Gene = gene_name,
        Taxonomy_group = group_name,
        Keyword = keyword,
        Matched_rows = sum(idx, na.rm = TRUE),
        t(summed),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
  )

  dplyr::bind_rows(aggregated_rows)
}


build_gene_sample_matrix <- function(combined_df, taxonomy_group) {

  metadata_cols <- c(
    "Gene",
    "Taxonomy_group",
    "Keyword",
    "Matched_rows"
  )

  sample_cols <- setdiff(colnames(combined_df), metadata_cols)

  sub_df <- combined_df %>%
    filter(Taxonomy_group == taxonomy_group) %>%
    select(Gene, all_of(sample_cols))

  if (nrow(sub_df) == 0) {
    stop("No data found for taxonomy group: ", taxonomy_group)
  }

  # If the same final gene name occurs more than once, sum abundance values.
  sub_df <- sub_df %>%
    group_by(Gene) %>%
    summarise(
      across(
        all_of(sample_cols),
        ~ sum(as.numeric(.x), na.rm = TRUE)
      ),
      .groups = "drop"
    )

  mat <- as.matrix(sub_df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- sub_df$Gene

  # Remove genes with insufficient variation because Spearman correlation
  # is undefined for constant vectors.
  keep <- apply(
    mat,
    1,
    function(x) {
      x <- x[is.finite(x)]
      length(x) >= 3 && stats::sd(x) > 0
    }
  )

  dropped <- rownames(mat)[!keep]
  mat <- mat[keep, , drop = FALSE]

  list(
    matrix = mat,
    dropped = dropped
  )
}


calculate_spearman_matrices <- function(mat) {

  genes <- rownames(mat)
  n_genes <- length(genes)

  rho_mat <- matrix(
    NA_real_,
    nrow = n_genes,
    ncol = n_genes,
    dimnames = list(genes, genes)
  )

  p_mat <- rho_mat
  n_mat <- rho_mat

  for (i in seq_len(n_genes)) {
    for (j in i:n_genes) {

      x <- mat[i, ]
      y <- mat[j, ]

      ok <- is.finite(x) & is.finite(y)
      n_pair <- sum(ok)

      n_mat[i, j] <- n_pair
      n_mat[j, i] <- n_pair

      if (i == j) {
        rho_mat[i, j] <- 1
        p_mat[i, j] <- 0
        next
      }

      if (
        n_pair < 3 ||
        stats::sd(x[ok]) == 0 ||
        stats::sd(y[ok]) == 0
      ) {
        next
      }

      test <- suppressWarnings(
        stats::cor.test(
          x[ok],
          y[ok],
          method = "spearman",
          exact = FALSE
        )
      )

      rho_value <- unname(test$estimate)
      p_value   <- test$p.value

      rho_mat[i, j] <- rho_value
      rho_mat[j, i] <- rho_value

      p_mat[i, j] <- p_value
      p_mat[j, i] <- p_value
    }
  }

  list(
    rho = rho_mat,
    p = p_mat,
    n = n_mat
  )
}


write_matrix_csv <- function(mat, file) {
  out <- data.frame(
    Gene = rownames(mat),
    mat,
    check.names = FALSE
  )

  utils::write.csv(
    out,
    file = file,
    row.names = FALSE,
    quote = TRUE
  )
}


build_complete_linkage_hclust <- function(cor_mat) {

  if (anyNA(cor_mat)) {
    stop(
      "The correlation matrix contains NA values and cannot be clustered. ",
      "Check whether one or more genes have insufficient overlapping or ",
      "variable abundance values."
    )
  }

  d <- stats::as.dist(1 - cor_mat)

  stats::hclust(
    d,
    method = "complete"
  )
}


make_correlation_heatmap <- function(
    cor_mat,
    title,
    show_legend = TRUE,
    fontsize = 6.0
) {

  hc <- build_complete_linkage_hclust(cor_mat)

  heat_colors <- grDevices::colorRampPalette(
    c("#2166AC", "#F7F7F7", "#B2182B")
  )(201)

  pheatmap::pheatmap(
    cor_mat,
    color = heat_colors,
    breaks = seq(-1, 1, length.out = 202),
    cluster_rows = hc,
    cluster_cols = hc,
    clustering_method = "complete",
    border_color = NA,
    legend = show_legend,
    main = title,
    fontsize = fontsize,
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    angle_col = 45,
    treeheight_row = 24,
    treeheight_col = 24,
    silent = TRUE
  )$gtable
}


save_grob <- function(
    grob,
    stem,
    width_cm,
    height_cm,
    dpi = 600
) {

  pdf_file  <- paste0(stem, ".pdf")
  png_file  <- paste0(stem, ".png")
  tiff_file <- paste0(stem, ".tiff")

  grDevices::pdf(
    pdf_file,
    width = width_cm / 2.54,
    height = height_cm / 2.54,
    useDingbats = FALSE
  )
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()

  grDevices::png(
    png_file,
    width = width_cm / 2.54,
    height = height_cm / 2.54,
    units = "in",
    res = dpi
  )
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()

  grDevices::tiff(
    tiff_file,
    width = width_cm / 2.54,
    height = height_cm / 2.54,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  grid::grid.newpage()
  grid::grid.draw(grob)
  grDevices::dev.off()
}


build_ut1_nr_data <- function(combined_df) {

  metadata_cols <- c(
    "Gene",
    "Taxonomy_group",
    "Keyword",
    "Matched_rows"
  )

  sample_cols <- setdiff(colnames(combined_df), metadata_cols)

  long_df <- combined_df %>%
    filter(
      Gene %in% c("UT1", "NR"),
      Taxonomy_group %in% scatter_taxa
    ) %>%
    select(
      Gene,
      Taxonomy_group,
      all_of(sample_cols)
    ) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "Abundance"
    ) %>%
    mutate(
      Abundance = suppressWarnings(as.numeric(Abundance))
    )

  wide_df <- long_df %>%
    select(
      Gene,
      Taxonomy_group,
      Sample,
      Abundance
    ) %>%
    pivot_wider(
      names_from = Gene,
      values_from = Abundance,
      values_fn = sum
    ) %>%
    filter(
      !is.na(UT1),
      !is.na(NR)
    ) %>%
    mutate(
      Taxonomy_group = factor(
        Taxonomy_group,
        levels = scatter_taxa
      ),
      log_NR = log10(NR + pseudocount),
      log_UT1 = log10(UT1 + pseudocount)
    )

  wide_df
}


calculate_ut1_nr_stats <- function(wide_df) {

  wide_df %>%
    group_by(Taxonomy_group) %>%
    group_modify(
      ~ {
        dat <- .x %>%
          filter(
            is.finite(NR),
            is.finite(UT1)
          )

        n_samples <- nrow(dat)

        if (
          n_samples < 3 ||
          sd(dat$NR) == 0 ||
          sd(dat$UT1) == 0
        ) {
          return(
            data.frame(
              n_samples = n_samples,
              rho = NA_real_,
              p_value = NA_real_,
              label = paste0("n = ", n_samples)
            )
          )
        }

        test <- suppressWarnings(
          cor.test(
            dat$NR,
            dat$UT1,
            method = "spearman",
            exact = FALSE
          )
        )

        rho_value <- unname(test$estimate)
        p_value   <- test$p.value

        p_label <- if (p_value < 0.001) {
          "P < 0.001"
        } else {
          paste0(
            "P = ",
            format.pval(
              p_value,
              digits = 2,
              eps = 0.001
            )
          )
        }

        data.frame(
          n_samples = n_samples,
          rho = rho_value,
          p_value = p_value,
          label = paste0(
            "\u03c1 = ",
            sprintf("%.2f", rho_value),
            "\n",
            p_label,
            "\nn = ",
            n_samples
          )
        )
      }
    ) %>%
    ungroup() %>%
    mutate(
      Taxonomy_group = factor(
        Taxonomy_group,
        levels = scatter_taxa
      )
    )
}


make_ut1_nr_scatter <- function(wide_df, stats_df) {

  ggplot(
    wide_df,
    aes(
      x = log_NR,
      y = log_UT1
    )
  ) +
    geom_point(
      size = 1.15,
      alpha = 0.70,
      color = "#34495E"
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.55,
      color = "#B2182B",
      fill = "#F4A6A6",
      alpha = 0.25
    ) +
    geom_text(
      data = stats_df,
      aes(
        x = -Inf,
        y = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = -0.06,
      vjust = 1.15,
      size = 2.5,
      color = "black"
    ) +
    facet_wrap(
      ~ Taxonomy_group,
      nrow = 1,
      scales = "free"
    ) +
    theme_classic(base_size = 8.5) +
    labs(
      title = expression(
        paste(
          italic("UT1"),
          "\u2013",
          italic("NR"),
          " Spearman correlation in MATOU_v1_metaT"
        )
      ),
      x = expression(
        log[10](
          italic(NR)~"abundance" + 1 %*% 10^-10
        )
      ),
      y = expression(
        log[10](
          italic(UT1)~"abundance" + 1 %*% 10^-10
        )
      )
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 9.0
      ),
      strip.background = element_rect(
        fill = "grey92",
        color = "grey70",
        linewidth = 0.30
      ),
      strip.text = element_text(
        face = "bold",
        size = 7.5
      ),
      axis.title = element_text(
        size = 7.0
      ),
      axis.text = element_text(
        size = 6.5
      ),
      panel.border = element_rect(
        color = "grey70",
        fill = NA,
        linewidth = 0.30
      ),
      panel.spacing = grid::unit(0.45, "lines"),
      plot.margin = margin(4, 4, 4, 4)
    )
}


# ==============================================================================
# 4. DISCOVER AND AGGREGATE ALL GENE-SPECIFIC ABUNDANCE MATRICES
# ==============================================================================

input_files <- list.files(
  raw_data_dir,
  pattern = input_file_pattern,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(input_files) == 0) {
  stop(
    "No .xlsx input files were found in: ",
    raw_data_dir
  )
}

input_files <- sort(input_files)

aggregated_list <- vector(
  mode = "list",
  length = length(input_files)
)

manifest_list <- vector(
  mode = "list",
  length = length(input_files)
)

for (i in seq_along(input_files)) {

  file <- input_files[i]
  gene_name <- infer_gene_name(file)

  aggregated <- aggregate_one_gene(
    file = file,
    gene_name = gene_name
  )

  aggregated_list[[i]] <- aggregated

  safe_gene <- safe_filename(gene_name)

  aggregated_file <- file.path(
    output_dir,
    "01_taxonomy_aggregated",
    paste0(
      safe_gene,
      "_taxonomy_merged_7levels.xlsx"
    )
  )

  openxlsx::write.xlsx(
    aggregated,
    aggregated_file,
    rowNames = FALSE,
    overwrite = TRUE
  )

  manifest_list[[i]] <- data.frame(
    Source_file = basename(file),
    Gene = gene_name,
    Number_of_input_rows = sum(aggregated$Matched_rows[
      aggregated$Taxonomy_group == "Eukaryota"
    ]),
    Number_of_samples = ncol(aggregated) - 4,
    stringsAsFactors = FALSE
  )
}

combined_df <- dplyr::bind_rows(aggregated_list)

# Order rows consistently by taxonomic level and then gene name.
combined_df <- combined_df %>%
  mutate(
    Taxonomy_group = factor(
      Taxonomy_group,
      levels = taxonomy_levels$Taxonomy_group
    )
  ) %>%
  arrange(
    Taxonomy_group,
    Gene
  ) %>%
  mutate(
    Taxonomy_group = as.character(Taxonomy_group)
  )

combined_file <- file.path(
  output_dir,
  "all_genes_taxonomy_merged_7levels.xlsx"
)

openxlsx::write.xlsx(
  combined_df,
  combined_file,
  rowNames = FALSE,
  overwrite = TRUE
)

manifest_df <- dplyr::bind_rows(manifest_list)

utils::write.csv(
  manifest_df,
  file.path(
    output_dir,
    "input_file_manifest.csv"
  ),
  row.names = FALSE
)

message(
  "Aggregated ",
  length(input_files),
  " gene-specific abundance matrices."
)


# ==============================================================================
# 5. PAIRWISE SPEARMAN CORRELATIONS FOR ALL SEVEN TAXONOMIC LEVELS
# ==============================================================================

correlation_results <- list()
heatmap_grobs <- list()
dropped_gene_records <- list()

for (taxon in taxonomy_levels$Taxonomy_group) {

  message("Calculating pairwise Spearman correlations: ", taxon)

  prepared <- build_gene_sample_matrix(
    combined_df,
    taxon
  )

  mat <- prepared$matrix

  if (nrow(mat) < 2) {
    stop(
      "Fewer than two variable genes were available for: ",
      taxon
    )
  }

  if (length(prepared$dropped) > 0) {
    dropped_gene_records[[taxon]] <- data.frame(
      Taxonomy_group = taxon,
      Gene = prepared$dropped,
      Reason = "Zero or insufficient abundance variation",
      stringsAsFactors = FALSE
    )
  }

  cor_res <- calculate_spearman_matrices(mat)

  correlation_results[[taxon]] <- cor_res

  taxon_safe <- safe_filename(taxon)

  write_matrix_csv(
    cor_res$rho,
    file.path(
      output_dir,
      "02_correlation_matrices",
      paste0(
        "Spearman_rho_",
        taxon_safe,
        ".csv"
      )
    )
  )

  write_matrix_csv(
    cor_res$p,
    file.path(
      output_dir,
      "02_correlation_matrices",
      paste0(
        "Spearman_pvalue_",
        taxon_safe,
        ".csv"
      )
    )
  )

  write_matrix_csv(
    cor_res$n,
    file.path(
      output_dir,
      "02_correlation_matrices",
      paste0(
        "Spearman_n_",
        taxon_safe,
        ".csv"
      )
    )
  )

  heatmap_grobs[[taxon]] <- make_correlation_heatmap(
    cor_mat = cor_res$rho,
    title = taxon,
    show_legend = TRUE,
    fontsize = 5.2
  )
}

if (length(dropped_gene_records) > 0) {

  dropped_df <- dplyr::bind_rows(
    dropped_gene_records
  )

  utils::write.csv(
    dropped_df,
    file.path(
      output_dir,
      "02_correlation_matrices",
      "genes_excluded_from_correlation_heatmaps.csv"
    ),
    row.names = FALSE
  )
}


# ==============================================================================
# 6. SUPPLEMENTARY FIG. S11a: SEVEN CORRELATION HEATMAPS
# ==============================================================================

# Layout:
#   Eukaryota          Stramenopiles       Bacillariophyta
#   Bacillariophycidae Coscinodiscophyceae Chaetocerotaceae
#   Thalassiosirales

s11a_layout <- matrix(
  c(
    1, 2, 3,
    4, 5, 6,
    7, NA, NA
  ),
  nrow = 3,
  byrow = TRUE
)

s11a_grob <- gridExtra::arrangeGrob(
  grobs = heatmap_grobs[
    taxonomy_levels$Taxonomy_group
  ],
  layout_matrix = s11a_layout,
  top = grid::textGrob(
    "Spearman correlations among nitrogen transporter and assimilation genes in MATOU_v1_metaT",
    gp = grid::gpar(
      fontsize = 11,
      fontface = "bold"
    )
  )
)

save_grob(
  grob = s11a_grob,
  stem = file.path(
    output_dir,
    "03_figures",
    "Supplementary_Fig_S11a_pairwise_Spearman_heatmaps"
  ),
  width_cm = 20,
  height_cm = 18
)


# ==============================================================================
# 7. FIG. 4a: THREE REPRESENTATIVE CORRELATION HEATMAPS
# ==============================================================================

fig4a_grobs <- lapply(
  seq_along(fig4a_taxa),
  function(i) {

    taxon <- fig4a_taxa[i]

    # To make Fig. 4a cleaner, show the color legend only on the final panel.
    make_correlation_heatmap(
      cor_mat = correlation_results[[taxon]]$rho,
      title = taxon,
      show_legend = (i == length(fig4a_taxa)),
      fontsize = 5.2
    )
  }
)

fig4a_grob <- gridExtra::arrangeGrob(
  grobs = fig4a_grobs,
  nrow = 1
)

save_grob(
  grob = fig4a_grob,
  stem = file.path(
    output_dir,
    "03_figures",
    "Fig4a_representative_Spearman_heatmaps"
  ),
  width_cm = 18,
  height_cm = 6.5
)


# ==============================================================================
# 8. SUPPLEMENTARY FIG. S11b: FOCUSED UT1-NR CORRELATION
# ==============================================================================

available_genes <- unique(combined_df$Gene)

if (!all(c("UT1", "NR") %in% available_genes)) {
  stop(
    "UT1 and/or NR were not found among the inferred gene names.\n",
    "Available genes include: ",
    paste(sort(available_genes), collapse = ", "),
    "\nUse gene_name_overrides if the input files use different names."
  )
}

ut1_nr_df <- build_ut1_nr_data(combined_df)

ut1_nr_stats <- calculate_ut1_nr_stats(
  ut1_nr_df
)

utils::write.csv(
  ut1_nr_stats %>%
    mutate(
      Taxonomy_group = as.character(
        Taxonomy_group
      )
    ) %>%
    select(
      Taxonomy_group,
      n_samples,
      rho,
      p_value
    ),
  file.path(
    output_dir,
    "02_correlation_matrices",
    "UT1_NR_Spearman_statistics.csv"
  ),
  row.names = FALSE
)

ut1_nr_plot <- make_ut1_nr_scatter(
  wide_df = ut1_nr_df,
  stats_df = ut1_nr_stats
)

ggsave(
  filename = file.path(
    output_dir,
    "03_figures",
    "Supplementary_Fig_S11b_UT1_NR_Spearman_scatter.pdf"
  ),
  plot = ut1_nr_plot,
  width = 12,
  height = 6,
  units = "cm",
  device = grDevices::cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "03_figures",
    "Supplementary_Fig_S11b_UT1_NR_Spearman_scatter.png"
  ),
  plot = ut1_nr_plot,
  width = 12,
  height = 6,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "03_figures",
    "Supplementary_Fig_S11b_UT1_NR_Spearman_scatter.tiff"
  ),
  plot = ut1_nr_plot,
  width = 12,
  height = 6,
  units = "cm",
  dpi = 600,
  compression = "lzw"
)


# ==============================================================================
# 9. OPTIONAL COMBINED SUPPLEMENTARY FIG. S11 LAYOUT
# ==============================================================================

# This reproduces the logical panel structure of Supplementary Fig. S11:
# six heatmaps in the first two rows, Thalassiosirales at lower left,
# and the UT1-NR scatter plot spanning the lower middle/right area.

s11_combined_layout <- matrix(
  c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 8
  ),
  nrow = 3,
  byrow = TRUE
)

s11_combined_grob <- gridExtra::arrangeGrob(
  grobs = c(
    heatmap_grobs[
      taxonomy_levels$Taxonomy_group
    ],
    list(ggplotGrob(ut1_nr_plot))
  ),
  layout_matrix = s11_combined_layout,
  top = grid::textGrob(
    "Co-variation of nitrogen transporter and assimilation homologues in Tara Oceans metatranscriptomes",
    gp = grid::gpar(
      fontsize = 11,
      fontface = "bold"
    )
  )
)

save_grob(
  grob = s11_combined_grob,
  stem = file.path(
    output_dir,
    "03_figures",
    "Supplementary_Fig_S11_complete_layout"
  ),
  width_cm = 20,
  height_cm = 18
)


# ==============================================================================
# 10. SAVE REPRODUCIBILITY INFORMATION
# ==============================================================================

session_file <- file.path(
  output_dir,
  "analysis_sessionInfo.txt"
)

capture.output(
  sessionInfo(),
  file = session_file
)


# ==============================================================================
# 11. FINAL SUMMARY
# ==============================================================================

message("")
message("==============================================================")
message("MATOU_v1_metaT analysis completed successfully.")
message("==============================================================")
message("Input gene files: ", length(input_files))
message("Combined abundance file: ", combined_file)
message("")
message("Main outputs:")
message("  - Per-gene seven-level abundance summaries")
message("  - Combined seven-level abundance matrix")
message("  - Spearman rho, exact P-value and n matrices for all seven taxa")
message("  - Fig. 4a three-taxonomy heatmaps")
message("  - Supplementary Fig. S11a seven-taxonomy heatmaps")
message("  - Supplementary Fig. S11b UT1-NR scatter plots")
message("  - Combined Supplementary Fig. S11 layout")
message("  - R session information")
message("")
message("Results directory: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
