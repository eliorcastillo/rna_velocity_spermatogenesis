# 04_velocity_dynamical.R
# Model 3: dynamical. Fits the full transcription/splicing/degradation
# kinetics per gene (alpha, beta, gamma) and a latent time per cell.
# No steady-state assumption. Slowest of the three.
# Output: rds/velo_dynamical.rds

source("00_setup.R")
sce_germ  <- readRDS("rds/sce_germ.rds")
hvgs_germ <- readRDS("rds/hvgs_germ.rds")

run_time <- system.time(
  velo_dynamical <- scvelo(
    sce_germ,
    subset.row    = hvgs_germ,
    assay.X       = "spliced",
    mode          = "dynamical", # <- the only line that differs between models
    scvelo.params = list(neighbors = list(n_neighbors = n_neighbors))
  )
)
run_time

stopifnot(identical(colnames(velo_dynamical), colnames(sce_germ)))

saveRDS(velo_dynamical, "rds/velo_dynamical.rds")
