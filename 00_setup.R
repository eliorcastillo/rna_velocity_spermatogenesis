# 00_setup.R
# Shared setup: packages, conflict preferences, output folder.
# Every other script starts with source("00_setup.R").

Sys.setenv(LANGUAGE = "en") # errors/warnings in English (en_EN.UTF-8 is not a valid locale on macOS)

# CRAN packages: course style, installed automatically if missing
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  conflicted, # detects and resolves function name conflicts
  tidyverse,  # data handling and ggplot2
  patchwork   # combining plots side by side
)

# Bioconductor packages: only checked, never installed from a script
pkgs_bioc <- c(
  "velociraptor", # R interface to scVelo
  "scRNAseq",     # Hermann spermatogenesis dataset
  "scuttle",      # normalisation
  "scran",        # highly variable genes
  "scater",       # PCA, t-SNE, plotting
  "org.Mm.eg.db"  # mouse gene annotation: Ensembl ID -> gene symbol
)
missing_pkgs <- pkgs_bioc[!vapply(pkgs_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing Bioconductor packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them with BiocManager::install(c(...)) and run again."
  )
}
pacman::p_load(char = pkgs_bioc)

# conflict_scout() # run interactively to check conflicts
conflicts_prefer(
  # tidyverse verbs used for tables
  dplyr::filter,
  dplyr::select, # org.Mm.eg.db loads AnnotationDbi, which has its own select()
  dplyr::rename,
  dplyr::first,
  dplyr::slice,
  dplyr::desc,
  dplyr::count,
  dplyr::combine,
  stats::lag,
  # set operations on gene lists
  base::intersect,
  base::union,
  base::setdiff,
  base::setequal,
  # matrix summaries on count matrices (sparse-aware)
  MatrixGenerics::rowMeans2,
  MatrixGenerics::rowVars,
  MatrixGenerics::rowMedians,
  .quiet = TRUE
)

# All saved objects go here
dir.create("rds", showWarnings = FALSE)

# Shared parameters, defined once so all three models use the same values
n_hvgs      <- 2000L
n_neighbors <- 30L
seed        <- 100L

# Stages: full names (as in the dataset, typo fixed) and short labels ------
stage_levels <- c(
  "Undifferentiated spermatogonia",
  "A1-A2 Differentiating spermatogonia",
  "A3-A4-In-B Differentiating spermatogonia",
  "Preleptotene spermatocytes",
  "Leptotene/Zygotene spermatocytes",
  "Pachytene spermatocytes",
  "Diplotene/Secondary spermatocytes",
  "Early Round spermatids",
  "Mid Round spermatids",
  "Late Round spermatids"
)
stage_short <- c(
  "Undiff SG", "A1-A2 SG", "A3-B SG", "Prelep", "Lep/Zyg",
  "Pach", "Dip/Sec", "Early RS", "Mid RS", "Late RS"
)
stage_map <- setNames(stage_short, stage_levels) # dictionary: long -> short

# Lineage palette: cold (spermatogonia) to warm (late spermatids)
stage_colours <- setNames(
  scales::viridis_pal(option = "turbo")(length(stage_short)),
  stage_short
)

# Helper: boxplot + individual cells per stage ------------------------------
# The number of cells per stage goes in the axis label ("Late RS (n=423)"),
# because early stages have very few cells.
plot_by_stage <- function(df, var, y_label) {
  counts <- df |> count(stage_short)
  axis_labels <- setNames(
    paste0(counts$stage_short, "\n(n=", counts$n, ")"),
    counts$stage_short
  )
  ggplot(df, aes(x = stage_short, y = .data[[var]], colour = stage_short)) +
    geom_boxplot(outlier.shape = NA, colour = "grey40") +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
    scale_x_discrete(labels = axis_labels) +
    scale_colour_manual(values = stage_colours, guide = "none") +
    labs(x = NULL, y = y_label) +
    theme_bw()
}

# Helper: t-SNE map coloured by any column of a per-cell table ---------------
plot_map <- function(df, var, title = var) {
  p <- ggplot(df, aes(x = tsne1, y = tsne2, colour = .data[[var]])) +
    geom_point(size = 0.8, alpha = 0.7) +
    labs(x = "t-SNE 1", y = "t-SNE 2", colour = NULL, title = title) +
    theme_bw()
  if (var == "stage_short") {
    p + scale_colour_manual(values = stage_colours) +
      guides(colour = guide_legend(override.aes = list(size = 3)))
  } else {
    p + scale_colour_viridis_c()
  }
}

# Helper: add velocity arrows (from gridVectors) to a map ------------------
add_arrows <- function(p, grid) {
  p + geom_segment(
    data = grid,
    aes(x = start.1, y = start.2, xend = end.1, yend = end.2),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.05, "inches")), linewidth = 0.3
  )
}

# Helper: add gene symbol and gene type to a per-gene table ------------------
# The dataset only has versioned Ensembl IDs (ENSMUSG00000038015.6).
# The version (.6) is removed before looking them up in org.Mm.eg.db.
# Genes not found in the database get NA.
add_gene_symbols <- function(df, id_col = "gene") {
  ens_ids <- sub("\\.\\d+$", "", df[[id_col]])
  lookup <- function(column) {
    suppressMessages( # hides the usual "1:many mapping" note; the first match is kept
      unname(AnnotationDbi::mapIds(
        org.Mm.eg.db,
        keys = ens_ids, keytype = "ENSEMBL", column = column, multiVals = "first"
      ))
    )
  }
  df |>
    mutate(
      symbol    = lookup("SYMBOL"),
      gene_type = lookup("GENETYPE"),
      .after = all_of(id_col)
    )
}