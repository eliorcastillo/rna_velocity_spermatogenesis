# 05_explore_steady_state.R
# Step 2, part 1: exploring the steady-state model.
#   Block 0: unspliced fraction per stage (raw counts, no model)
#   Block 1: the funnel (cells and genes)
#   Block 2: per-gene results (gamma, fit quality)
#   Block 3: per-cell results (confidence, velocity length)
#   Block 4: the order (velocity pseudotime per stage)
# Reads rds files only: fast, can be rerun as often as needed.
# Output: figures/ and rds/cells_steady.rds, rds/genes_steady.rds

source("00_setup.R")
sce_germ    <- readRDS("rds/sce_germ.rds")
hvgs_germ   <- readRDS("rds/hvgs_germ.rds")
velo_steady <- readRDS("rds/velo_steady.rds")
stopifnot(identical(colnames(velo_steady), colnames(sce_germ)))

dir.create("figures", showWarnings = FALSE)

# Plotting helpers (stage_short, stage_map, stage_colours, plot_by_stage)
# now live in 00_setup.R, so every script can use them.


# Block 0: unspliced fraction per stage (no model) --------------------------
# Hypothesis: if transcription shuts down in late spermatids, their
# fraction of unspliced (nascent) RNA should be the lowest.

qc_s     <- perCellQCMetrics(sce_germ, assay.type = "spliced")
qc_u     <- perCellQCMetrics(sce_germ, assay.type = "unspliced")
qc_s_hvg <- perCellQCMetrics(sce_germ[hvgs_germ, ], assay.type = "spliced")
qc_u_hvg <- perCellQCMetrics(sce_germ[hvgs_germ, ], assay.type = "unspliced")

cells <- tibble(
  cell        = colnames(sce_germ),
  stage       = sce_germ$stage,
  stage_short = factor(stage_map[as.character(sce_germ$stage)], levels = stage_short),
  spliced     = qc_s$sum,
  unspliced   = qc_u$sum,
  total       = spliced + unspliced,
  frac_unspliced     = unspliced / total,                                   # all genes
  frac_unspliced_hvg = qc_u_hvg$sum / (qc_s_hvg$sum + qc_u_hvg$sum),        # the 2,000 HVGs
  tsne1 = reducedDim(sce_germ, "TSNE")[, 1],
  tsne2 = reducedDim(sce_germ, "TSNE")[, 2]
)

# Table: one row per stage
cells |>
  group_by(stage_short) |>
  summarise(
    n_cells                   = n(),
    median_total_counts       = median(total),
    median_frac_unspliced     = round(median(frac_unspliced), 3),
    median_frac_unspliced_hvg = round(median(frac_unspliced_hvg), 3)
  ) |>
  print(n = Inf)

p0a <- plot_by_stage(cells, "frac_unspliced", "Unspliced fraction (all genes)")
p0b <- plot_by_stage(cells, "total", "Total counts per cell (spliced + unspliced)") +
  scale_y_log10()
p0a / p0b
ggsave("figures/00_unspliced_fraction.pdf", p0a / p0b, width = 9, height = 9)


# Block 1: the funnel --------------------------------------------------------

funnel_steady <- tibble(
  step = c("Cells", "Genes in dataset", "HVGs given to the model", "Velocity genes"),
  n = c(
    ncol(velo_steady),
    nrow(sce_germ),
    nrow(velo_steady),
    sum(rowData(velo_steady)$velocity_genes, na.rm = TRUE)
  )
)
funnel_steady


# Block 2: per-gene results --------------------------------------------------

# What did scVelo store per gene, and what do gene names look like?
colnames(rowData(velo_steady))
head(rownames(velo_steady))

genes_steady <- as.data.frame(rowData(velo_steady)) |>
  as_tibble(rownames = "gene")

# Gamma: the steady-state ratio u/s for each gene
summary(genes_steady$velocity_gamma)

p2 <- genes_steady |>
  filter(velocity_gamma > 0) |>
  ggplot(aes(x = log10(velocity_gamma), fill = velocity_genes)) +
  geom_histogram(bins = 50) +
  labs(x = "log10(gamma)", y = "Genes", fill = "Velocity gene") +
  theme_bw()
p2
ggsave("figures/02_steady_gamma.pdf", p2, width = 7, height = 4)

# Gene symbols and types (function defined in 00_setup.R)
genes_steady <- add_gene_symbols(genes_steady)
sum(is.na(genes_steady$symbol))                       # genes without a symbol

# Are velocity genes enriched in some gene type (e.g. lncRNAs)?
genes_steady |> count(velocity_genes, gene_type) |> print(n = Inf)

# Atypical genes: highest gamma (more unspliced than expected)
genes_steady |>
  arrange(desc(velocity_gamma)) |>
  select(symbol, gene_type, velocity_gamma, velocity_r2, velocity_genes) |>
  head(10)

# Best-fitting velocity genes (highest R2 of the steady-state line)
genes_steady |>
  filter(velocity_genes) |>
  arrange(desc(velocity_r2)) |>
  select(symbol, gene_type, velocity_gamma, velocity_r2) |>
  head(20)


# Block 3: per-cell results --------------------------------------------------

cells <- cells |>
  mutate(
    pseudotime = velo_steady$velocity_pseudotime,
    confidence = velo_steady$velocity_confidence,
    length     = velo_steady$velocity_length,
    root       = velo_steady$root_cells,
    end        = velo_steady$end_points
  )

cells |>
  group_by(stage_short) |>
  summarise(
    n_cells           = n(),
    median_confidence = round(median(confidence), 3),
    median_length     = round(median(length), 3)
  ) |>
  print(n = Inf)

p3a <- plot_by_stage(cells, "confidence", "Velocity confidence")
p3b <- plot_by_stage(cells, "length", "Velocity length")
p3a / p3b
ggsave("figures/03_steady_confidence_length.pdf", p3a / p3b, width = 9, height = 9)


# Block 4: the order ---------------------------------------------------------

order_steady <- cells |>
  group_by(stage_short) |>
  summarise(
    n_cells           = n(),
    median_pseudotime = round(median(pseudotime), 3)
  )
print(order_steady, n = Inf)

# Does the median pseudotime increase with the known lineage order?
# FALSE means at least one stage is out of order.
!is.unsorted(order_steady$median_pseudotime)

p4 <- plot_by_stage(cells, "pseudotime", "Velocity pseudotime (steady-state)")
p4
ggsave("figures/04_steady_pseudotime.pdf", p4, width = 9, height = 5)


# Save for the comparison step -----------------------------------------------
saveRDS(cells,        "rds/cells_steady.rds")
saveRDS(genes_steady, "rds/genes_steady.rds")