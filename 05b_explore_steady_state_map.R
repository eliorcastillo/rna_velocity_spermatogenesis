# 05b_explore_steady_state_map.R
# Step 2, part 2: the steady-state model on the map.
#   Block 5: root and end cells (where does the model think it starts/ends?)
#   Block 6: velocity arrows on the t-SNE
#   Block 7: the isolated island of mid round spermatids
# Needs: rds/cells_steady.rds (from 05_explore_steady_state.R)
# Output: figures/ and rds/grid_steady.rds, rds/cells_steady.rds (updated)

source("00_setup.R")
sce_germ    <- readRDS("rds/sce_germ.rds")
velo_steady <- readRDS("rds/velo_steady.rds")
cells       <- readRDS("rds/cells_steady.rds")
stopifnot(identical(cells$cell, colnames(velo_steady)))


# Block 5: root and end cells -------------------------------------------------
# root/end are probabilities (0-1). A cell counts as root/end if > 0.5.

roots_ends <- cells |>
  group_by(stage_short) |>
  summarise(
    n_cells   = n(),
    n_root    = sum(root > 0.5),
    prop_root = round(mean(root > 0.5), 2),
    n_end     = sum(end > 0.5),
    prop_end  = round(mean(end > 0.5), 2)
  )
print(roots_ends, n = Inf)

p5 <- plot_map(cells, "root", "Root probability") |
  plot_map(cells, "end", "End-point probability")
p5
ggsave("figures/05_steady_roots_ends.pdf", p5, width = 12, height = 5)


# Block 6: velocity arrows ----------------------------------------------------
# embedVelocity() projects each cell's velocity onto the t-SNE (runs Python);
# gridVectors() averages them on a grid so the map stays readable.

emb_steady  <- embedVelocity(reducedDim(sce_germ, "TSNE"), velo_steady)
grid_steady <- gridVectors(sce_germ, emb_steady, use.dimred = "TSNE")
saveRDS(grid_steady, "rds/grid_steady.rds") # reused in the comparison step

p6a <- add_arrows(plot_map(cells, "stage_short", "Stages"), grid_steady)
p6b <- add_arrows(plot_map(cells, "pseudotime", "Velocity pseudotime (steady-state)"), grid_steady)
p6a | p6b
ggsave("figures/06_steady_arrows.pdf", p6a | p6b, width = 14, height = 6)


# Block 7: the island of mid round spermatids -------------------------------
# Selected by position on the map. The box is read from the t-SNE plot:
# CHECK the highlighted map below and adjust the limits if needed.

island_box <- list(x = c(-10, 5), y = c(-13, -2))

cells <- cells |>
  mutate(
    island = stage_short == "Mid RS" &
      between(tsne1, island_box$x[1], island_box$x[2]) &
      between(tsne2, island_box$y[1], island_box$y[2])
  )
sum(cells$island)

# Visual check: is the selection exactly the island?
p7a <- ggplot(cells, aes(tsne1, tsne2, colour = island)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey80", `TRUE` = "red")) +
  labs(title = "Island selection") +
  theme_bw()
p7a

# Island vs the rest of mid round spermatids
cells |>
  filter(stage_short == "Mid RS") |>
  group_by(island) |>
  summarise(
    n_cells           = n(),
    median_total      = median(total),
    median_frac_unspl = round(median(frac_unspliced), 3),
    median_confidence = round(median(confidence), 3),
    median_length     = round(median(length), 2),
    median_pseudotime = round(median(pseudotime), 3),
    prop_root         = round(mean(root > 0.5), 2),
    prop_end          = round(mean(end > 0.5), 2)
  )

ggsave("figures/07_island_selection.pdf", p7a, width = 7, height = 6)

# Save the table again, now with the island column
saveRDS(cells, "rds/cells_steady.rds")









cells |>
  filter(stage_short == "Mid RS") |>
  group_by(island) |>
  summarise(prop_root = round(mean(root > 0.5), 2), prop_end = round(mean(end > 0.5), 2))
