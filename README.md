# RNA velocity in mouse spermatogenesis

Can RNA velocity recover the known order of spermatogenesis from a single snapshot, and where does it fail?

Spermatogenesis is a continuous, mostly linear differentiation with a well-known order of stages, so the correct answer is known in advance. This makes it a good test case for learning what RNA velocity can and cannot tell us.

## Data

- **Experiment:** Hermann et al. 2018, *Cell Reports*. Mouse spermatogenic cells, 10x Chromium 3′ scRNA-seq.
- **Spliced/unspliced counts:** requantified by Soneson et al. 2021, *PLoS Computational Biology*.
- **Access:** `scRNAseq::HermannSpermatogenesisData()` (Bioconductor). Data are downloaded by the package; nothing is stored in this repo.

Germ cells only (Sertoli cells and unannotated cells removed): 1,876 cells in 10 stages, from undifferentiated spermatogonia to late round spermatids. Early stages have very few cells (1–10 per stage).

## Scripts

Run from the repository root, in order. Every script starts with `source("00_setup.R")`.

| Script | What it does | Status |
|---|---|---|
| `00_setup.R` | Packages, shared parameters, stage labels, plotting helpers | ✅ |
| `01_prepare_data.R` | Germ cells, normalisation, 2,000 HVGs, PCA and t-SNE | ✅ |
| `02_velocity_steady_state.R` | scVelo, steady-state model | ✅ |
| `03_velocity_stochastic.R` | scVelo, stochastic model | ⚠️ not run (Python incompatibility) |
| `04_velocity_dynamical.R` | scVelo, dynamical model | ✅ fitted; exploration pending |
| `05_explore_steady_state.R` | Unspliced fraction, genes, cells, pseudotime order | ✅ |
| `05b_explore_steady_state_map.R` | Root/end cells, velocity arrows, mid round spermatid island | ✅ |

Intermediate objects are saved to `rds/` (not tracked) and figures to `figures/`.

## Key findings (steady-state model)

1. **Unspliced fraction drops in late round spermatids**, before any model is fitted (`figures/00_unspliced_fraction.pdf`). This is the expected signature of transcriptional shutdown: unspliced RNA is cleared quickly, while stabilised spliced mRNA persists.
2. **The last stage is inverted.** Root cells are placed among late round spermatids and their velocity pseudotime is low, so the model reads the end of the lineage as its start (`figures/04_steady_pseudotime.pdf`, `figures/05_steady_roots_ends.pdf`).
3. **Interpretation:** velocity correctly detects that genes are being switched off in late round spermatids, but the true next state (elongating spermatids) is not in the dataset. The closest "future-like" cells are therefore earlier stages.
4. **An isolated island of mid round spermatids** on the t-SNE concentrates end-point probability (`figures/06_steady_arrows.pdf`, `figures/07_island_selection.pdf`). Whether it is biology or a technical subgroup is still open.

## Caveats

- Preprocessing choices (counting tool, spliced/unspliced assignment) change velocity results, as shown by Soneson et al. 2021 on this same dataset.
- t-SNE distorts distances; arrows on the map are a summary, not a measurement.
- Early stages are represented by very few cells.
- This is a learning exercise, not a new result.

## Next steps

- Explore the dynamical model (latent time, per-gene kinetics) and compare with steady-state.
- Resolve the Python issue and run the stochastic model.

## References

- La Manno et al. 2018, *Nature*. RNA velocity of single cells.
- Bergen et al. 2020, *Nature Biotechnology*. scVelo.
- Hermann et al. 2018, *Cell Reports*. The mammalian spermatogenesis single-cell transcriptome.
- Soneson et al. 2021, *PLoS Computational Biology*. Preprocessing choices affect RNA velocity results for droplet scRNA-seq data.
