# 03_velocity_stochastic.R
# Model 2: stochastic. Same steady-state idea, but the fit also uses
# second-order moments (variances/covariances of u and s), not only means.
# Output: rds/velo_stochastic.rds

# source("00_setup.R")
# sce_germ  <- readRDS("rds/sce_germ.rds")
# hvgs_germ <- readRDS("rds/hvgs_germ.rds")
# 
# run_time <- system.time(
#   velo_stochastic <- scvelo(
#     sce_germ,
#     subset.row    = hvgs_germ,
#     assay.X       = "spliced",
#     mode          = "stochastic", # <- the only line that differs between models
#     scvelo.params = list(neighbors = list(n_neighbors = n_neighbors))
#   )
# )
# run_time
# 
# stopifnot(identical(colnames(velo_stochastic), colnames(sce_germ)))
# 
# saveRDS(velo_stochastic, "rds/velo_stochastic.rds")
# 
# 
# reticulate::py_last_error()                         # 1. the full Python traceback
# reticulate::import("numpy")$`__version__`           # 2. numpy version
# reticulate::import("scvelo")$`__version__`          # 3. scVelo version


#Python incompatibilites