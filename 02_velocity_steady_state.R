# 02_velocity_steady_state.R
# Model 1: steady-state (deterministic). One gamma per gene, fitted on
# the extreme quantiles of the phase portrait. scVelo's default.
# Output: rds/velo_steady.rds

source("00_setup.R")
sce_germ  <- readRDS("rds/sce_germ.rds")
hvgs_germ <- readRDS("rds/hvgs_germ.rds")

run_time <- system.time(
  velo_steady <- scvelo(
    sce_germ,
    subset.row    = hvgs_germ,
    assay.X       = "spliced",
    mode          = "steady_state", # <- the only line that differs between models
    scvelo.params = list(neighbors = list(n_neighbors = n_neighbors))
  )
)
run_time

colData(velo_steady)
sum(rowData(velo_steady)$velocity_genes)

# Cells must be in the same order as sce_germ, otherwise every comparison is wrong
stopifnot(identical(colnames(velo_steady), colnames(sce_germ)))

saveRDS(velo_steady, "rds/velo_steady.rds")
