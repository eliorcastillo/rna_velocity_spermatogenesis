# 01_prepare_data.R
# Builds the single input object used by all three velocity models:
# germ cells only, normalised, 2,000 HVGs, one PCA and one t-SNE.
# Output: rds/sce_germ.rds, rds/hvgs_germ.rds, rds/dec_germ.rds

source("00_setup.R")

# 1. Load the dataset (spliced + unspliced counts, Hermann et al. 2018) ----
sce_full <- HermannSpermatogenesisData()
sce_full
table(sce_full$celltype, useNA = "ifany")

# 2. Keep the germ-cell lineage only ----
# Sertoli cells act as a second, unrelated "root" and NA cells have no label.
keep <- !is.na(sce_full$celltype) & sce_full$celltype != "Sertoli cells"

sce_germ <- sce_full[, keep]
sce_germ
colData(sce_full)

# 3. Stage as an ordered factor (lineage order), fixing the "DIplotene" typo ----
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
sce_germ$stage <- factor(
  sub("^DIplotene", "Diplotene", sce_germ$celltype),
  levels = stage_levels
)
stopifnot(!anyNA(sce_germ$stage)) # every cell must match a known stage
table(sce_germ$stage)
ncol(sce_germ) # The total number of cells

colData(sce_germ) #Here we can check the new "stage" column.

class(sce_germ$celltype)   # "character"
class(sce_germ$stage)      # "factor"

# 4. Normalise the spliced counts ----
# The unspliced counts are normalised later, inside scvelo().
sce_germ <- logNormCounts(sce_germ, assay.type = "spliced")

# 5. Highly variable genes, computed on germ cells only ----
dec_germ  <- modelGeneVar(sce_germ)
hvgs_germ <- getTopHVGs(dec_germ, n = n_hvgs)
length(hvgs_germ)

# 6. One shared map for all models ----
set.seed(seed)
sce_germ <- runPCA(sce_germ, subset_row = hvgs_germ)
sce_germ <- runTSNE(sce_germ, dimred = "PCA", perplexity = 30)
plotTSNE(sce_germ, colour_by = "stage")

# 7. Save ----
saveRDS(sce_germ,  "rds/sce_germ.rds")
saveRDS(hvgs_germ, "rds/hvgs_germ.rds")
saveRDS(dec_germ,  "rds/dec_germ.rds")
