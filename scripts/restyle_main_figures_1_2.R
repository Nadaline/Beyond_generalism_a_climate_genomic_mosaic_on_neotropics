
# 1. Load packages, saved results and plotting functions ----

library(ggplot2)
library(patchwork)
library(vegan)
library(maps)
library(terra)

project <- normalizePath(getwd())
results <- Sys.getenv("AMEIVA_RESULTS_DIR", file.path(project, "results", "RDA-2027-14134"))
figures <- file.path(project, "results", "figures", "RDA-2027-14134")
manuscript_figures <- file.path(project, "manuscript_figures_final")
environment_file <- Sys.getenv("AMEIVA_ENV_FILE", file.path("data", "local", "Ameiva_Variables_rda.csv"))

purple_dark <- "#32134F"
purple <- "#5B2A86"
lavender <- "#B79AD2"
orange_light <- "#F7C58D"
orange <- "#C76100"
cream <- "#FFF7EF"
snmf_colours <- c(`1` = purple_dark, `2` = lavender, `3` = orange_light, `4` = orange)
candidate_colours <- c(`1` = purple, `2` = "#F28E2B")

theme_ameiva <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey91", linewidth = .3),
      plot.title = element_text(face = "bold", colour = purple_dark, size = base_size + 2, margin = margin(b = 5)),
      axis.title = element_text(colour = "grey25"),
      legend.title = element_text(face = "bold", colour = purple_dark),
      legend.position = "bottom"
    )
}

save_plate <- function(x, name, width, height) {
  for (directory in c(figures, manuscript_figures)) {
    ggsave(file.path(directory, paste0(name, ".png")), x, width = width, height = height, dpi = 600, bg = "white")
    ggsave(file.path(directory, paste0(name, ".tiff")), x, width = width, height = height, dpi = 600,
           compression = "lzw", bg = "white")
  }
}

# 2. Plot sampling geography, genome-wide PCA and sNMF coefficients ----

q <- read.csv(file.path(results, "03_sNMF_Q.csv"), stringsAsFactors = FALSE)
pca <- read.csv(file.path(results, "04_PCA_14134_scores.csv"), stringsAsFactors = FALSE)
pca_model <- readRDS(file.path(results, "04_PCA_14134.rds"))
variables <- read.csv(environment_file, check.names = FALSE, stringsAsFactors = FALSE)
q <- q[match(pca$InternalID, q$InternalID), ]
variables <- variables[match(pca$InternalID, variables$InternalID), ]
stopifnot(identical(q$InternalID, pca$InternalID), identical(variables$InternalID, pca$InternalID), nrow(q) == 245L)

# Calculate percentages from the full PCA eigenspectrum, not the exported score columns.
pca_axis_fraction <- pca_model$CA$eig / sum(pca_model$CA$eig)

q_matrix <- as.matrix(q[, c("V1", "V2", "V3", "V4")])
q_matrix <- sweep(q_matrix, 1, rowSums(q_matrix), "/")
pca$snmf_cluster <- factor(max.col(q_matrix), levels = 1:4)
structure <- data.frame(InternalID = pca$InternalID, longitude = variables$long, latitude = variables$lat,
                        snmf_cluster = pca$snmf_cluster)
world <- map_data("world")

p1_map <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = cream, colour = "grey76", linewidth = .2) +
  geom_point(data = structure, aes(longitude, latitude, fill = snmf_cluster), shape = 21, colour = "white", stroke = .35, size = 2.65,
             position = position_jitter(width = .10, height = .10, seed = 2027)) +
  scale_fill_manual(values = snmf_colours, guide = "none") +
  coord_quickmap(xlim = c(-82, -32), ylim = c(-35, 12)) +
  labs(title = "Continental sampling and ancestry", x = "Longitude", y = "Latitude") +
  theme_ameiva() + theme(panel.grid = element_blank())

p1_pca <- ggplot(pca, aes(PC1, PC2, colour = snmf_cluster)) +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "grey90", linewidth = .3) +
  geom_point(size = 2.35, alpha = .9) +
  stat_ellipse(aes(group = snmf_cluster), linewidth = .45, level = .80, show.legend = FALSE) +
  scale_colour_manual(values = snmf_colours, guide = "none") +
  labs(title = "Genomic structure across Ameiva ameiva",
       x = paste0("PC1 (", round(100 * pca_axis_fraction[1], 1), "% variance)"),
       y = paste0("PC2 (", round(100 * pca_axis_fraction[2], 1), "% variance)")) +
  theme_ameiva()

structure_order <- order(max.col(q_matrix), variables$long)
q_long <- do.call(rbind, lapply(1:4, function(k) data.frame(
  sample_order = seq_along(structure_order),
  sNMF_cluster = factor(k, levels = 1:4),
  ancestry = q_matrix[structure_order, k]
)))
p1_q <- ggplot(q_long, aes(sample_order, ancestry, fill = sNMF_cluster)) +
  geom_col(width = 1.01, colour = NA) +
  scale_fill_manual(values = snmf_colours, name = "sNMF cluster") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(title = "Continuous ancestry at sNMF K = 4", x = "Individuals ordered by sNMF cluster and longitude", y = "Ancestry proportion") +
  theme_ameiva() + theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")

figure_1 <- (p1_map + p1_pca) / p1_q +
  plot_layout(heights = c(1, .67)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13, colour = purple_dark))
save_plate(figure_1, "Main_Figure_1_structure", 13, 10)
file.copy(file.path(manuscript_figures, "Main_Figure_1_structure.png"), file.path(manuscript_figures, "Figure_1_structure.png"), overwrite = TRUE)

# 3. Plot partial RDA, candidate-genotype clusters and projected axes ----

# Extract site and biplot scores from the saved partial RDA.
rda_geo <- readRDS(file.path(results, "07_RDA_PCA_geo.rds"))
rda_sites <- as.data.frame(scores(rda_geo, display = "sites", choices = 1:2, scaling = 2))
rda_arrows <- as.data.frame(scores(rda_geo, display = "bp", choices = 1:2, scaling = 2))
names(rda_sites)[1:2] <- c("RDA1", "RDA2")
names(rda_arrows)[1:2] <- c("RDA1", "RDA2")
rda_sites$snmf_cluster <- pca$snmf_cluster
rda_arrows$environment <- rownames(rda_arrows)
arrow_multiplier <- .70 * min(max(abs(rda_sites$RDA1)) / max(abs(rda_arrows$RDA1)),
                                max(abs(rda_sites$RDA2)) / max(abs(rda_arrows$RDA2)))
rda_arrows$RDA1 <- rda_arrows$RDA1 * arrow_multiplier
rda_arrows$RDA2 <- rda_arrows$RDA2 * arrow_multiplier
rda_axis_fraction <- rda_geo$CCA$eig / sum(rda_geo$CCA$eig)

p2_rda <- ggplot(rda_sites, aes(RDA1, RDA2, colour = snmf_cluster)) +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = .3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = .3) +
  geom_point(size = 2.05, alpha = .84) +
  geom_segment(data = rda_arrows, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), inherit.aes = FALSE,
               arrow = grid::arrow(length = grid::unit(2.0, "mm")), colour = orange, linewidth = .65) +
  geom_text(data = rda_arrows, aes(RDA1, RDA2, label = environment), inherit.aes = FALSE,
            colour = orange, fontface = "bold", size = 3) +
  scale_colour_manual(values = snmf_colours, guide = "none") + coord_equal() +
  labs(title = "Structure-controlled association",
       x = paste0("RDA1 (", round(100 * rda_axis_fraction[1], 1), "% constrained)"),
       y = paste0("RDA2 (", round(100 * rda_axis_fraction[2], 1), "% constrained)")) +
  theme_ameiva()

candidate_file <- file.path(results, "09_LFMM2_RDA_candidate_genotype_PCA", "03_RDA_LFMM2_shared_candidate_genotype_PCA_scores.csv")
candidate <- read.csv(candidate_file, stringsAsFactors = FALSE)
candidate$genotype_cluster <- factor(candidate$genotype_cluster, levels = 1:2)

p2_pca <- ggplot(candidate, aes(PC1, PC2, colour = genotype_cluster)) +
  geom_hline(yintercept = 0, colour = "grey90", linewidth = .3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = .3) +
  geom_point(size = 2.25, alpha = .9) + scale_colour_manual(values = candidate_colours, name = "Genotype cluster") +
  labs(title = "Linear candidates genotype cluster", x = "Candidate-genotype PC1", y = "Candidate-genotype PC2") +
  theme_ameiva() + theme(plot.margin = margin(l = 18, r = 5.5, b = 5.5, t = 5.5))

p2_map <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = cream, colour = "grey76", linewidth = .2) +
  geom_point(data = candidate, aes(longitude, latitude, fill = genotype_cluster), shape = 21, colour = "white", stroke = .35, size = 2.8,
             position = position_jitter(width = .12, height = .12, seed = 4027)) +
  scale_fill_manual(values = candidate_colours, guide = "none") +
  coord_quickmap(xlim = c(-82, -32), ylim = c(-35, 12)) +
  labs(title = "Spatial distribution of candidate-genotype clusters", x = "Longitude", y = "Latitude") +
  theme_ameiva() + theme(panel.grid = element_blank())

# Existing present-day RDA axes use diverging colour scales based on
# symmetric absolute 2nd–98th percentile limits from each full-resolution axis.
rda_landscape_file <- file.path(results, "spatial_offset", "18_RDA_present_axes.tif")
rda_landscape <- terra::rast(rda_landscape_file)
rda_display_limits <- lapply(1:2, function(i) {
  axis_quantiles <- quantile(terra::values(rda_landscape[[i]]), probs = c(.02, .98), na.rm = TRUE)
  max(abs(axis_quantiles)) * c(-1, 1)
})
rda1_display_limit <- rda_display_limits[[1]]
rda2_display_limit <- rda_display_limits[[2]]
rda_landscape_plot <- terra::aggregate(rda_landscape, fact = 4, fun = mean, na.rm = TRUE)
rda_landscape_df <- as.data.frame(rda_landscape_plot[[1]], xy = TRUE, na.rm = TRUE)
names(rda_landscape_df)[3] <- "RDA1"
p2_landscape <- ggplot(rda_landscape_df, aes(x, y, fill = RDA1)) +
  geom_raster(interpolate = FALSE) +
  geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE,
               fill = NA, colour = "grey72", linewidth = .18) +
  coord_quickmap(xlim = c(-82, -32), ylim = c(-35, 12), expand = FALSE) +
  scale_fill_gradient2(low = purple_dark, mid = cream, high = orange,
                       midpoint = 0, limits = rda1_display_limit, oob = scales::squish,
                       breaks = sort(c(rda1_display_limit, 0)),
                       labels = function(x) formatC(x, format = "f", digits = 1),
                       name = "RDA1 score") +
  labs(title = "Linear Climate-associated Genomic Landscape", x = "Longitude", y = "Latitude") +
  theme_ameiva() + theme(panel.grid = element_blank())

rda2_landscape_df <- as.data.frame(rda_landscape_plot[[2]], xy = TRUE, na.rm = TRUE)
names(rda2_landscape_df)[3] <- "RDA2"
p2_landscape_rda2 <- ggplot(rda2_landscape_df, aes(x, y, fill = RDA2)) +
  geom_raster(interpolate = FALSE) +
  geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE,
               fill = NA, colour = "grey72", linewidth = .18) +
  coord_quickmap(xlim = c(-82, -32), ylim = c(-35, 12), expand = FALSE) +
  scale_fill_gradient2(low = purple_dark, mid = cream, high = orange,
                       midpoint = 0, limits = rda2_display_limit, oob = scales::squish,
                       breaks = sort(c(rda2_display_limit, 0)),
                       labels = function(x) formatC(x, format = "f", digits = 1),
                       name = "RDA2 score") +
  labs(title = "Linear Climate-associated Genomic Landscape", x = "Longitude", y = "Latitude") +
  theme_ameiva() + theme(panel.grid = element_blank())

figure_2_top <- p2_rda + p2_pca + plot_layout(ncol = 2, widths = c(1, 1))
figure_2_bottom <- p2_map + p2_landscape + p2_landscape_rda2 +
  plot_layout(ncol = 3, widths = c(1, 1, 1))
figure_2 <- figure_2_top / figure_2_bottom +
  plot_layout(heights = c(.85, 1.35)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13, colour = purple_dark))
save_plate(figure_2, "Main_Figure_2_linear_adaptive_mosaic", 16, 15)
file.copy(file.path(manuscript_figures, "Main_Figure_2_linear_adaptive_mosaic.png"), file.path(manuscript_figures, "Figure_2_linear_adaptive_mosaic.png"), overwrite = TRUE)

message("Figures 1 and 2 restyled | saved analysis objects only | no analysis rerun")
