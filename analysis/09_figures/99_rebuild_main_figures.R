###############################################################################
# Reconstrucao leve das quatro pranchas principais e dos graficos finais de
# offset. Este script NAO refaz as analises: le somente resultados auditados do
# HD externo e imagens ggplot produzidas por RDA-2027.R.
###############################################################################

library(ggplot2)
library(cowplot)
library(png)
library(terra)
library(maps)

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"
hd <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation"
results <- file.path(hd, "results", "RDA-2027-14134")
figures <- file.path(project, "results", "figures", "RDA-2027-14134")
spatial <- file.path(results, "spatial_offset")
rf_spatial <- file.path(spatial, "rdaforest_original")

stopifnot(dir.exists(results), dir.exists(figures))

purple_dark <- "#2F1647"
purple <- "#5B2A86"
purple_light <- "#B99AD0"
orange_dark <- "#C76100"
orange <- "#F28E2B"
orange_light <- "#F8C58F"
cream <- "#FFF7EF"
grey <- "#BDBDBD"

sequential <- c(
  cream, orange_light, orange, "#A95572", purple, purple_dark
)

theme_ameiva <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey91", linewidth = .3),
      plot.title = element_text(face = "bold", colour = purple_dark),
      plot.subtitle = element_text(colour = "grey35"),
      strip.text = element_text(face = "bold", colour = purple_dark),
      strip.background = element_rect(fill = "#F3EAF8", colour = NA),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(colour = "grey20")
    )
}

idv_false_colour <- function(red, green, blue) {
  channels <- cbind(red, green, blue)
  channels[!is.finite(channels)] <- 0
  channels[channels < 0] <- 0
  channels[channels > 1] <- 1
  anchors <- t(col2rgb(c(purple, orange, cream)) / 255)
  weights <- channels + .08
  mixed <- weights %*% anchors / rowSums(weights)
  rgb(mixed[, 1], mixed[, 2], mixed[, 3])
}

image_panel <- function(filename) {
  path <- file.path(figures, filename)
  stopifnot(file.exists(path))
  image_grob <- grid::rasterGrob(
    readPNG(path),
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    interpolate = TRUE
  )
  ggdraw() + draw_grob(image_grob, x = 0, y = 0, width = 1, height = 1)
}

save_plate <- function(plot, filename, width, height) {
  ggsave(
    file.path(figures, paste0(filename, ".png")),
    plot, width = width, height = height,
    dpi = 400, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(figures, paste0(filename, ".tiff")),
    plot, width = width, height = height,
    dpi = 400, bg = "white", compression = "lzw", limitsize = FALSE
  )
}

# Graficos 26 e 28 sao refeitos diretamente das tabelas espaciais salvas.
offsets <- read.csv(file.path(
  rf_spatial, "21_RDAforest_original_offset_pixels.csv"
))
offsets$period <- factor(
  offsets$period, levels = c("2041-2070", "2071-2100")
)
plot_offset <- ggplot(offsets, aes(x, y, colour = offset)) +
  geom_point(shape = 15, size = .62) +
  facet_wrap(~period, nrow = 1) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
  scale_colour_gradientn(colours = sequential) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "RDAforest genomic displacement",
    subtitle = "Original gen_offset_oj; relative present-day q90 units",
    x = "Longitude", y = "Latitude", colour = "Relative\ndisplacement"
  )
ggsave(
  file.path(figures, "26_RDAforest_original_genomic_offset.png"),
  plot_offset, width = 11, height = 5.5, dpi = 400, bg = "white"
)

domain <- read.csv(file.path(
  rf_spatial, "21_RDAforest_projection_domain_pixels.csv"
))
domain$period <- factor(
  domain$period, levels = c("Present", "2041-2070", "2071-2100")
)
domain$status <- factor(
  domain$status, levels = c("inside", "clamped", "excluded")
)
plot_domain <- ggplot(
  subset(domain, period != "Present"),
  aes(x, y, colour = status)
) +
  geom_point(shape = 15, size = .55) +
  facet_wrap(~period, nrow = 1) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
  scale_colour_manual(values = c(
    inside = purple, clamped = orange, excluded = grey
  )) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "RDAforest projection-domain audit",
    subtitle = "Purple: training range; orange: clamped; grey: excluded",
    x = "Longitude", y = "Latitude", colour = "Projection status"
  )
ggsave(
  file.path(figures, "28_RDAforest_projection_domain_audit.png"),
  plot_domain, width = 11, height = 5.5, dpi = 400, bg = "white"
)

# Reconstroi a comparacao 27 pela extracao dos rasters GF e RDA exatamente nas
# coordenadas dos pixels RDAforest.
compare_period <- function(period_name, gf_file, rda_file) {
  period_data <- subset(offsets, period == period_name)
  gf <- rast(file.path(spatial, gf_file))
  rda <- rast(file.path(spatial, rda_file))
  points <- vect(period_data, geom = c("x", "y"), crs = crs(gf))
  extracted <- extract(c(gf, rda), points)
  result <- data.frame(
    x = period_data$x,
    y = period_data$y,
    RDAforest = period_data$offset,
    Gradient_Forest = extracted[[2]],
    RDA = extracted[[3]],
    period = period_name
  )
  result[complete.cases(result), , drop = FALSE]
}

offset_comparison <- rbind(
  compare_period(
    "2041-2070", "19_GF_offset_mid_century.tif",
    "19_RDA_offset_mid_century.tif"
  ),
  compare_period(
    "2071-2100", "19_GF_offset_late_century.tif",
    "19_RDA_offset_late_century.tif"
  )
)
offset_rho <- tapply(
  seq_len(nrow(offset_comparison)), offset_comparison$period,
  function(index) cor(
    offset_comparison$RDAforest[index],
    offset_comparison$Gradient_Forest[index],
    method = "spearman"
  )
)
plot_comparison <- ggplot(
  offset_comparison, aes(RDAforest, Gradient_Forest)
) +
  geom_point(alpha = .10, size = .45, colour = purple) +
  stat_summary_bin(
    fun = median, geom = "line", bins = 35,
    colour = orange, linewidth = 1
  ) +
  facet_wrap(~period, scales = "free") +
  theme_ameiva() +
  labs(
    title = "Agreement between genomic-displacement indices",
    subtitle = sprintf(
      "Spearman rho: %.2f (2041-2070) and %.2f (2071-2100)",
      offset_rho[["2041-2070"]], offset_rho[["2071-2100"]]
    ),
    x = "RDAforest displacement (present q90 units)",
    y = "Gradient Forest transformed-space displacement"
  )
ggsave(
  file.path(figures, "27_RDAforest_GF_offset_comparison.png"),
  plot_comparison, width = 10, height = 5.5, dpi = 400, bg = "white"
)

# Reconstroi somente o mapa presente do RDAforest a partir do ajuste combinado.
grid_data <- readRDS(file.path(rf_spatial, "21_RDAforest_spatial_grid.rds"))
combined_fit <- readRDS(file.path(
  rf_spatial, "21_RDAforest_combined_ordinationJackknife.rds"
))
period_n <- nrow(grid_data$xy)
retained_global <- which(combined_fit$goodrows)
present_prediction_rows <- which(retained_global <= period_n)
present_good <- combined_fit$goodrows[seq_len(period_n)]
present_predictions <- combined_fit$predictions.direct[
  present_prediction_rows, , drop = FALSE
]
present_xy <- grid_data$xy[present_good, , drop = FALSE]
stopifnot(nrow(present_predictions) == nrow(present_xy))

prediction_pca <- prcomp(
  present_predictions, center = TRUE, scale. = TRUE
)
present_scores <- prediction_pca$x[, 1:3, drop = FALSE]
present_limits <- apply(
  present_scores, 2, quantile, probs = c(.02, .98), na.rm = TRUE
)
present_scores01 <- sweep(present_scores, 2, present_limits[1, ], "-")
present_scores01 <- sweep(
  present_scores01, 2, present_limits[2, ] - present_limits[1, ], "/"
)
present_scores01[present_scores01 < 0] <- 0
present_scores01[present_scores01 > 1] <- 1
present_xy$brand_colour <- idv_false_colour(
  present_scores01[, 1], present_scores01[, 2], present_scores01[, 3]
)

plot_rda_forest_present <- ggplot(
  present_xy, aes(x, y, colour = brand_colour)
) +
  annotation_borders("world", fill = cream, colour = "grey78") +
  geom_point(shape = 15, size = .65) +
  scale_colour_identity() +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Present adaptive landscape — RDAforest",
    subtitle = "Genome-wide residual covariance; project false colour",
    x = "Longitude", y = "Latitude"
  )

# Titulo abreviado para manter o painel candidato legivel na prancha tripla.
mosaic <- read.csv(file.path(
  results, "13_adaptive_mosaic_genotype_clusters.csv"
))
mosaic$genotype_cluster <- factor(mosaic$genotype_cluster)
plot_candidate_pca <- ggplot(
  mosaic, aes(PC1, PC2, colour = genotype_cluster)
) +
  geom_point(size = 2.7, alpha = .9) +
  scale_colour_manual(values = c(purple, orange)) +
  theme_ameiva() +
  labs(
    title = "Candidate-genotype mosaic",
    subtitle = "PCA of 104 SNPs supported jointly by RDA and LFMM",
    x = "PC1", y = "PC2", colour = "Genotype cluster"
  )
ggsave(
  file.path(figures, "05_adaptive_mosaic_genotype_PCA.png"),
  plot_candidate_pca, width = 7, height = 5.5,
  dpi = 400, bg = "white"
)

rank01 <- function(x) (rank(x, ties.method = "average") - 1) / (length(x) - 1)
offset_comparison$linear_rank <- ave(
  offset_comparison$RDA, offset_comparison$period, FUN = rank01
)
offset_comparison$GF_rank <- ave(
  offset_comparison$Gradient_Forest, offset_comparison$period, FUN = rank01
)
offset_comparison$RDAforest_rank <- ave(
  offset_comparison$RDAforest, offset_comparison$period, FUN = rank01
)
offset_comparison$nonlinear_rank <- rowMeans(
  offset_comparison[, c("GF_rank", "RDAforest_rank")]
)
offset_comparison$agreement <- 1 - abs(
  offset_comparison$linear_rank - offset_comparison$nonlinear_rank
)
offset_comparison$interpretation <- "concordant intermediate"
offset_comparison$interpretation[
  offset_comparison$linear_rank >= .67 & offset_comparison$nonlinear_rank >= .67
] <- "concordant high"
offset_comparison$interpretation[
  offset_comparison$linear_rank <= .33 & offset_comparison$nonlinear_rank <= .33
] <- "concordant low"
offset_comparison$interpretation[
  offset_comparison$linear_rank - offset_comparison$nonlinear_rank > .25
] <- "linear higher"
offset_comparison$interpretation[
  offset_comparison$nonlinear_rank - offset_comparison$linear_rank > .25
] <- "nonlinear higher"
offset_comparison$interpretation <- factor(
  offset_comparison$interpretation,
  levels = c(
    "concordant low", "concordant intermediate", "concordant high",
    "linear higher", "nonlinear higher"
  )
)
concordance_late <- subset(offset_comparison, period == "2071-2100")
write.csv(
  offset_comparison,
  file.path(figures, "35_linear_nonlinear_offset_concordance_source.csv"),
  row.names = FALSE
)

plot_linear_rank <- ggplot(
  concordance_late, aes(x, y, colour = linear_rank)
) +
  geom_point(shape = 15, size = .62) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
  scale_colour_gradientn(colours = sequential, limits = c(0, 1)) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Linear climate association",
    subtitle = "Partial-RDA late-century offset rank",
    x = "Longitude", y = "Latitude", colour = "Continental\nrank"
  )

plot_GF_rank <- ggplot(
  concordance_late, aes(x, y, colour = GF_rank)
) +
  geom_point(shape = 15, size = .62) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
  scale_colour_gradientn(colours = sequential, limits = c(0, 1)) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Candidate nonlinear association",
    subtitle = "Gradient Forest late-century offset rank",
    x = "Longitude", y = "Latitude", colour = "Continental\nrank"
  )

plot_RDAforest_rank <- ggplot(
  concordance_late, aes(x, y, colour = RDAforest_rank)
) +
  geom_point(shape = 15, size = .62) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
  scale_colour_gradientn(colours = sequential, limits = c(0, 1)) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Genome-wide nonlinear association",
    subtitle = "RDAforest late-century offset rank",
    x = "Longitude", y = "Latitude", colour = "Continental\nrank"
  )

concordance_colours <- c(
  "concordant low" = purple_dark,
  "concordant intermediate" = purple_light,
  "concordant high" = orange_dark,
  "linear higher" = purple,
  "nonlinear higher" = orange
)
plot_concordance <- ggplot(
  concordance_late, aes(x, y, colour = interpretation)
) +
  geom_point(shape = 15, size = .62) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
  scale_colour_manual(values = concordance_colours, drop = FALSE) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Linear–nonlinear concordance",
    subtitle = "Agreement of relative late-century displacement",
    x = "Longitude", y = "Latitude", colour = "Spatial reading"
  )

# Prancha 1: dois paineis superiores e uma barra sNMF ocupando toda a base.
figure_1_top <- plot_grid(
  image_panel("32_population_sampling_map.png"),
  image_panel("02_PCA_sNMF_structure.png"),
  labels = c("a", "b"), label_size = 16, ncol = 2
)
figure_1 <- plot_grid(
  figure_1_top,
  image_panel("33_sNMF_ancestry_barplot.png"),
  labels = c("", "c"), label_size = 16,
  ncol = 1, rel_heights = c(1, .66)
)

figure_2 <- plot_grid(
  image_panel("31_partial_RDA_genotype_response.png"),
  image_panel("05_adaptive_mosaic_genotype_PCA.png"),
  image_panel("06_adaptive_mosaic_map.png"),
  labels = c("a", "b", "c"), label_size = 16,
  ncol = 3, rel_widths = c(1.12, 1, 1.12)
)

figure_3 <- plot_grid(
  image_panel("03_Gradient_Forest_importance.png"),
  image_panel("18_RDAforest_importance_uncertainty.png"),
  image_panel("34_candidate_functional_context.png"),
  labels = c("a", "b", "c"), label_size = 16, ncol = 3
)

figure_4 <- plot_grid(
  plot_linear_rank,
  plot_GF_rank,
  plot_RDAforest_rank,
  plot_concordance,
  labels = c("a", "b", "c", "d"), label_size = 16,
  ncol = 2
)

save_plate(figure_1, "Main_Figure_1_structure", 13, 10)
save_plate(figure_2, "Main_Figure_2_linear_adaptive_mosaic", 16, 6.2)
save_plate(figure_3, "Main_Figure_3_polygenic_functional_architecture", 16, 6.2)
save_plate(figure_4, "Main_Figure_4_landscape_and_offset", 14, 10)

write.csv(
  data.frame(
    figure = paste0("Main_Figure_", 1:4),
    panels = c(3, 3, 3, 4),
    biological_role = c(
      "continental sampling, genomic PCA and continuous K=4 ancestry",
      "partial RDA with allelic response and the RDA-and-LFMM candidate mosaic",
      "candidate GF, genome-wide RDAforest and functional context",
      "late-century RDA, GF and RDAforest ranks plus linear-nonlinear concordance"
    ),
    png = c(
      "Main_Figure_1_structure.png",
      "Main_Figure_2_linear_adaptive_mosaic.png",
      "Main_Figure_3_polygenic_functional_architecture.png",
      "Main_Figure_4_landscape_and_offset.png"
    ),
    tiff = c(
      "Main_Figure_1_structure.tiff",
      "Main_Figure_2_linear_adaptive_mosaic.tiff",
      "Main_Figure_3_polygenic_functional_architecture.tiff",
      "Main_Figure_4_landscape_and_offset.tiff"
    )
  ),
  file.path(figures, "MAIN_FIGURE_REGISTRY.csv"),
  row.names = FALSE
)

cat("Four main plates, offset diagnostics and concordance rebuilt successfully.\n")
