###############################################################################
# PRESENT-DAY NONLINEAR GENOMIC LANDSCAPES — POST-PROCESSING ONLY
# Reuses completed controlled GF blocks and final saved RDAforest predictions.
# No forest, jackknife, candidate screen, fold, or RDA is refitted.
###############################################################################

#### 1. PACKAGES, PATHS AND CURRENT OBJECT LINEAGE ####

library(terra)
library(ggplot2)
library(patchwork)
library(maps)
library(gradientForest)
terraOptions(memfrac = .20, tempdir = file.path(tempdir(), "ameiva_nonlinear_landscapes"))

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"
results <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation/results/RDA-2027-14134"
out <- file.path(project, "results", "nonlinear_present_landscapes")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

gf_dir <- file.path(results, "gradient_forest_all_14134", "controlled_PC1_PC3_lon_lat")
gf_blocks <- sort(list.files(gf_dir, "^GF_controlled_block_.*\\.rds$", full.names = TRUE))
rda_template_file <- file.path(results, "spatial_offset", "18_RDA_present_axes.tif")
rdaforest_dir <- file.path(results, "spatial_offset", "rdaforest_original")
rdaforest_grid_file <- file.path(rdaforest_dir, "21_RDAforest_spatial_grid.rds")
rdaforest_spatial_file <- file.path(rdaforest_dir, "21_RDAforest_combined_ordinationJackknife.rds")
rdaforest_signature_file <- file.path(rdaforest_dir, "21_RDAforest_spatial_signature.csv")
rdaforest_final_file <- file.path(results, "12_RDAforest_fit_v2610_nonquarter_geoonly_ecdf_v2_pcs8.rds")
environment_file <- "/home/nadaline/Documents/Phd/RADseq/Ameiva_Variables_rda.csv"
climate_dir <- "/projects/rstudio-export/Ameiva"
stopifnot(length(gf_blocks) == 32L, all(file.exists(c(rda_template_file, rdaforest_grid_file,
  rdaforest_spatial_file, rdaforest_signature_file, rdaforest_final_file, environment_file))))

predictors <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio12", "bio13", "bio15")
variables <- read.csv(environment_file, check.names = FALSE, stringsAsFactors = FALSE)
world <- map_data("world")
template <- rast(rda_template_file)[[1]]

lineage <- c(
  "CONTROLLED GF: 32 saved GF_controlled_block_*_of_32.rds objects under gradient_forest_all_14134/controlled_PC1_PC3_lon_lat/.",
  "These are the current 14,134-SNP controlled forests (PC1-PC3 + longitude + latitude removed before GF).",
  paste0("RDAFOREST: ", rdaforest_final_file, "; spatial predictions from ", rdaforest_spatial_file, "."),
  paste0("RDAforest grid signature: ", rdaforest_signature_file, "."),
  "Deprecated candidate-only GF products and non-final RDAforest directories were not used."
)
writeLines(lineage, file.path(out, "nonlinear_present_landscape_object_audit.txt"))

#### 2. CONTROLLED GF TRANSFORMATION OF PRESENT-DAY CLIMATE ####

# Each block has the same ten ECDF predictors and a different SNP subset. The
# mean of its already fitted cumulative-importance curves is a post-processing
# summary of the whole controlled 14,134-SNP GF panel; no model is fit here.
curve_grid <- seq(0, 1, length.out = 1001)
gf_curve_summary <- lapply(predictors, function(p) {
  y <- vapply(gf_blocks, function(f) {
    block <- readRDS(f)
    curve <- gradientForest::cumimp(block, p)
    approx(curve$x, curve$y, xout = curve_grid, rule = 2, ties = "ordered")$y
  }, numeric(length(curve_grid)))
  data.frame(predictor = p, ecdf = curve_grid, mean_cumimp = rowMeans(y),
             sd_cumimp = apply(y, 1, sd))
})
gf_curve_summary <- do.call(rbind, gf_curve_summary)
write.csv(gf_curve_summary, file.path(out, "GF_controlled_mean_cumulative_curves.csv"), row.names = FALSE)

# The raw CHELSA layers are only read and aligned to the pre-existing RDA
# present-day template/mask. ECDFs are calculated from the original 245
# environmental observations, matching the controlled GF input convention.
gf_layers <- lapply(predictors, function(p) {
  message("Transforming controlled GF predictor: ", p)
  raw <- rast(file.path(climate_dir, paste0("CHELSA_", p, "_1981-2010_SA.tif")))
  raw <- crop(raw, ext(-85, -30, -35, 15))
  raw <- aggregate(raw, fact = 5, fun = mean, na.rm = TRUE)
  stopifnot(compareGeom(raw, template, stopOnError = FALSE))
  raw <- mask(raw, template)
  observed <- variables[[paste0("X.CHELSA_", p, "_1981.2010_SA")]]
  observed <- sort(observed[is.finite(observed)])
  empirical_probability <- seq_along(observed) / length(observed)
  curve <- gf_curve_summary[gf_curve_summary$predictor == p, ]
  ecdf_raster <- app(raw, function(v) approx(observed, empirical_probability,
                                              xout = v, rule = 2, ties = "ordered")$y)
  transformed <- app(ecdf_raster, function(v) approx(curve$ecdf, curve$mean_cumimp,
                                                      xout = v, rule = 2)$y)
  names(transformed) <- p
  transformed
})
gf_transformed <- rast(gf_layers)
names(gf_transformed) <- predictors

#### 3. VISUALIZATION-ONLY GF PCA ####

set.seed(2027)
gf_sample <- spatSample(gf_transformed, size = 100000, method = "random", na.rm = TRUE,
                        as.df = TRUE, values = TRUE)
gf_sample <- gf_sample[complete.cases(gf_sample[, predictors]), predictors, drop = FALSE]
gf_pca <- prcomp(gf_sample, center = TRUE, scale. = TRUE)
saveRDS(gf_pca, file.path(out, "GF_present_visualization_PCA.rds"))
gf_variance <- 100 * gf_pca$sdev^2 / sum(gf_pca$sdev^2)

gf_standardized <- lapply(predictors, function(p) (gf_transformed[[p]] - gf_pca$center[p]) / gf_pca$scale[p])
gf_pc1 <- sum(rast(gf_standardized) * gf_pca$rotation[predictors, 1])
gf_pc2 <- sum(rast(gf_standardized) * gf_pca$rotation[predictors, 2])
names(gf_pc1) <- "GF_landscape_score"
names(gf_pc2) <- "GF_landscape_PC2"
writeRaster(gf_pc1, file.path(out, "GF_present_landscape.tif"), overwrite = TRUE)
writeRaster(gf_pc2, file.path(out, "GF_present_landscape_PC2.tif"), overwrite = TRUE)

#### 4. SAVED PRESENT-DAY RDAFOREST PREDICTIONS ####

signature <- read.csv(rdaforest_signature_file, check.names = FALSE)
stopifnot(signature$selected_variables == "bio13+bio3+bio2", signature$top_pcs == 8)
grid <- readRDS(rdaforest_grid_file)
combined <- readRDS(rdaforest_spatial_file)
final_fit <- readRDS(rdaforest_final_file)
period_n <- nrow(grid$xy)
present_good <- combined$goodrows[seq_len(period_n)]
present_direct <- combined$predictions.direct[seq_len(sum(present_good)), , drop = FALSE]
stopifnot(ncol(present_direct) == 8L, nrow(present_direct) == sum(present_good), final_fit$failed.reps == 0)

# The saved spatial workflow used an unsaved visualization PCA of direct
# predictions. Recreating that PCA from exactly the same saved prediction
# matrix is permitted visualization-only post-processing.
rdaforest_pca <- prcomp(present_direct, center = TRUE, scale. = TRUE)
saveRDS(rdaforest_pca, file.path(out, "RDAforest_present_visualization_PCA.rds"))
rdaforest_variance <- 100 * rdaforest_pca$sdev^2 / sum(rdaforest_pca$sdev^2)
rdaforest_scores <- predict(rdaforest_pca, present_direct)[, 1:2, drop = FALSE]
rdaforest_xy <- grid$xy[present_good, , drop = FALSE]
rdaforest_df <- data.frame(x = rdaforest_xy$x, y = rdaforest_xy$y,
                            RDAforest_landscape_score = rdaforest_scores[, 1],
                            RDAforest_landscape_PC2 = rdaforest_scores[, 2])
rdaforest_pc1 <- rast(rdaforest_df[, c("x", "y", "RDAforest_landscape_score")], type = "xyz",
                       crs = crs(template))
rdaforest_pc2 <- rast(rdaforest_df[, c("x", "y", "RDAforest_landscape_PC2")], type = "xyz",
                       crs = crs(template))
writeRaster(rdaforest_pc1, file.path(out, "RDAforest_present_landscape.tif"), overwrite = TRUE)
writeRaster(rdaforest_pc2, file.path(out, "RDAforest_present_landscape_PC2.tif"), overwrite = TRUE)

#### 5. TAIL AUDITS: SAME P10/P90 AND ONE-DEGREE LOGIC AS RDA ####

tail_audit <- function(x, name) {
  v <- values(x)
  p10 <- unname(quantile(v, .10, na.rm = TRUE))
  p90 <- unname(quantile(v, .90, na.rm = TRUE))
  low <- ifel(x <= p10, 1, 0)
  high <- ifel(x >= p90, 1, 0)
  names(low) <- paste0(name, "_tail_P10")
  names(high) <- paste0(name, "_tail_P90")
  list(
    summary = data.frame(method = name, P10 = p10, P90 = p90,
                         finite_cells = sum(is.finite(v)), lower_tail_cells = sum(v <= p10, na.rm = TRUE),
                         upper_tail_cells = sum(v >= p90, na.rm = TRUE)),
    low = low, high = high
  )
}

gf_tail <- tail_audit(gf_pc1, "Gradient Forest")
rdaforest_tail <- tail_audit(rdaforest_pc1, "RDAforest")
writeRaster(gf_tail$low, file.path(out, "GF_tail_P10.tif"), overwrite = TRUE)
writeRaster(gf_tail$high, file.path(out, "GF_tail_P90.tif"), overwrite = TRUE)
writeRaster(rdaforest_tail$low, file.path(out, "RDAforest_tail_P10.tif"), overwrite = TRUE)
writeRaster(rdaforest_tail$high, file.path(out, "RDAforest_tail_P90.tif"), overwrite = TRUE)
write.csv(gf_tail$summary, file.path(out, "GF_tail_audit.csv"), row.names = FALSE)
write.csv(rdaforest_tail$summary, file.path(out, "RDAforest_tail_audit.csv"), row.names = FALSE)

tail_grid <- function(low, high, method) {
  # Existing RDA grid is 0.0416667 degrees; factor 24 gives the same 1-degree comparison grid.
  fact <- max(1L, round(1 / res(low)[1]))
  make <- function(x, label) {
    y <- aggregate(x, fact = fact, fun = mean, na.rm = TRUE)
    d <- as.data.frame(y, xy = TRUE, na.rm = TRUE)
    names(d)[3] <- "tail_fraction"
    d$method <- method; d$tail <- label; d
  }
  rbind(make(low, "Lower tail (≤ P10)"), make(high, "Upper tail (≥ P90)"))
}
gf_tail_grid <- tail_grid(gf_tail$low, gf_tail$high, "Gradient Forest")
# RDAforest predictions are already approximately 0.5-degree grid cells;
# aggregate by two for the same one-degree tail-display resolution.
rda_fact <- max(1L, round(1 / res(rdaforest_tail$low)[1]))
make_rda_grid <- function(x, label) {
  y <- aggregate(x, fact = rda_fact, fun = mean, na.rm = TRUE)
  d <- as.data.frame(y, xy = TRUE, na.rm = TRUE); names(d)[3] <- "tail_fraction"
  d$method <- "RDAforest"; d$tail <- label; d
}
rdaforest_tail_grid <- rbind(make_rda_grid(rdaforest_tail$low, "Lower tail (≤ P10)"),
                             make_rda_grid(rdaforest_tail$high, "Upper tail (≥ P90)"))
all_tail_grid <- rbind(gf_tail_grid, rdaforest_tail_grid)
write.csv(all_tail_grid, file.path(out, "nonlinear_present_tail_concentration_1deg.csv"), row.names = FALSE)

tail_map <- ggplot(all_tail_grid, aes(x, y, fill = tail_fraction)) +
  geom_raster() + geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE,
                                fill = NA, colour = "grey68", linewidth = .18) +
  coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15), expand = FALSE) +
  facet_grid(method ~ tail) + scale_fill_gradient(low = "white", high = "#5B2A86", name = "Fraction of\n1° cell in tail") +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "bottom",
                                       strip.text = element_text(face = "bold", colour = "#32134F")) +
  labs(x = "Longitude", y = "Latitude", title = "Geographic concentration of nonlinear landscape tails")
ggsave(file.path(out, "Figure_nonlinear_present_tail_audit.png"), tail_map, width = 11, height = 7, dpi = 300, bg = "white")

#### 6. SECONDARY AXES AND INTERNAL SPATIAL CORRELATIONS ####

gf_pc2_tail <- tail_audit(gf_pc2, "Gradient Forest PC2")
rdaforest_pc2_tail <- tail_audit(rdaforest_pc2, "RDAforest PC2")
write.csv(rbind(gf_pc2_tail$summary, rdaforest_pc2_tail$summary),
          file.path(out, "nonlinear_secondary_axis_tail_audit.csv"), row.names = FALSE)

rda_axes <- rast(rda_template_file)
common_xy <- rdaforest_xy
gf_at_grid <- extract(gf_pc1, common_xy, ID = FALSE)[, 1]
rda1_at_grid <- extract(rda_axes[[1]], common_xy, ID = FALSE)[, 1]
rda2_at_grid <- extract(rda_axes[[2]], common_xy, ID = FALSE)[, 1]
correlation_data <- data.frame(GF = gf_at_grid, RDAforest = rdaforest_scores[, 1], RDA1 = rda1_at_grid, RDA2 = rda2_at_grid)
correlation_data <- correlation_data[complete.cases(correlation_data), ]
cmat <- cor(correlation_data, method = "spearman")
correlation_out <- data.frame()
for (i in seq_len(ncol(cmat) - 1)) for (j in (i + 1):ncol(cmat)) correlation_out <- rbind(correlation_out,
  data.frame(axis_x = colnames(cmat)[i], axis_y = colnames(cmat)[j], signed_spearman_rho = cmat[i, j],
             absolute_spearman_rho = abs(cmat[i, j]), common_cells = nrow(correlation_data)))
write.csv(correlation_out, file.path(out, "nonlinear_present_landscape_spatial_correlations.csv"), row.names = FALSE)

#### 7. TWO-PANEL MANUSCRIPT-STYLE MAP ####

display_map <- function(x, title, legend_label) {
  q <- quantile(values(x), c(.02, .98), na.rm = TRUE)
  limits <- max(abs(q)) * c(-1, 1)
  xd <- aggregate(x, fact = 4, fun = mean, na.rm = TRUE)
  d <- as.data.frame(xd, xy = TRUE, na.rm = TRUE); names(d)[3] <- "score"
  ggplot(d, aes(x, y, fill = score)) + geom_raster() +
    geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE,
                 fill = NA, colour = "grey72", linewidth = .18) +
    geom_point(data = variables, aes(long, lat), inherit.aes = FALSE,
               shape = 21, fill = "white", colour = "#32134F", stroke = .25, size = .8, alpha = .7) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15), expand = FALSE) +
    scale_fill_gradient2(low = "#32134F", mid = "#FFF7EF", high = "#F28E2B", midpoint = 0,
                         limits = limits, oob = scales::squish, name = legend_label,
                         breaks = sort(c(limits, 0)), labels = function(z) formatC(z, format = "f", digits = 1)) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 10) + theme(panel.grid = element_blank(),
      plot.title = element_text(face = "bold", colour = "#32134F"),
      legend.position = "bottom", legend.title = element_text(face = "bold", colour = "#32134F"))
}

p_gf <- display_map(gf_pc1, "Gradient Forest", "GF landscape score") + labs(tag = "a")
p_rdaforest <- display_map(rdaforest_pc1, "RDAforest", "RDAforest genomic score") + labs(tag = "b")
figure <- p_gf + p_rdaforest + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13, colour = "#32134F"))
ggsave(file.path(out, "Figure_present_GF_RDAforest_v1.pdf"), figure, width = 12, height = 6.5, device = cairo_pdf)
ggsave(file.path(out, "Figure_present_GF_RDAforest_v1.png"), figure, width = 12, height = 6.5, dpi = 600, bg = "white")

#### 8. NOTES FOR MANUAL REGIONAL INTERPRETATION ####

audit <- c(
  "PRESENT-DAY NONLINEAR LANDSCAPE AUDIT", "",
  lineage, "",
  "GRADIENT FOREST",
  "Displayed quantity: visualization-only PC1 of the existing controlled-GF cumulative-transform surface.",
  sprintf("Variance represented: PC1 %.2f%%; PC2 %.2f%%; PC1+PC2 %.2f%%.", gf_variance[1], gf_variance[2], sum(gf_variance[1:2])),
  sprintf("P10: %.6f; P90: %.6f; finite cells: %d.", gf_tail$summary$P10, gf_tail$summary$P90, gf_tail$summary$finite_cells),
  "Sign orientation: retained as returned by the visualization PCA; not flipped.",
  "", "RDAFOREST",
  "Displayed quantity: visualization-only PC1 of existing final RDAforest present-day direct predictions.",
  sprintf("Variance represented: PC1 %.2f%%; PC2 %.2f%%; PC1+PC2 %.2f%%.", rdaforest_variance[1], rdaforest_variance[2], sum(rdaforest_variance[1:2])),
  sprintf("P10: %.6f; P90: %.6f; finite cells: %d.", rdaforest_tail$summary$P10, rdaforest_tail$summary$P90, rdaforest_tail$summary$finite_cells),
  "Sign orientation: retained as returned by the visualization PCA; not flipped.",
  "", "Secondary-axis note: P10/P90 numerical audits are in nonlinear_secondary_axis_tail_audit.csv; map inspection is required before naming regions.",
  "Internal correlations only: nonlinear_present_landscape_spatial_correlations.csv.",
  "No model was refit and no underlying predicted value was modified."
)
writeLines(audit, file.path(out, "nonlinear_present_landscape_audit.txt"))

cat(paste(c(
  "GF object: 32 controlled whole-panel GF block objects.",
  "RDAforest object: final nonquarter_geoonly_ecdf_v2_pcs8 fit plus saved rdaforest_original present predictions.",
  "Visualization PCA: required for both methods; no signs flipped.",
  sprintf("GF PC1 variance %.2f%% | RDAforest PC1 variance %.2f%%.", gf_variance[1], rdaforest_variance[1]),
  sprintf("GF P10/P90 %.6f / %.6f | RDAforest P10/P90 %.6f / %.6f.", gf_tail$summary$P10, gf_tail$summary$P90, rdaforest_tail$summary$P10, rdaforest_tail$summary$P90),
  paste0("Final maps: ", file.path(out, "Figure_present_GF_RDAforest_v1.pdf"), " ; ", file.path(out, "Figure_present_GF_RDAforest_v1.png"))
), collapse = "\n"), "\n")
