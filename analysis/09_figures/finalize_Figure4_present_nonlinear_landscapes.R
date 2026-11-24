###############################################################################
# FIGURE 4 — PRESENT-DAY NONLINEAR LANDSCAPES
# Display-only revision: existing rasters and raw visualization-PC orientations.
###############################################################################

#### 1. INPUTS AND LINEAGE ####

library(terra)
library(ggplot2)
library(patchwork)
library(maps)
library(mgcv)

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"
out <- file.path(project, "results", "nonlinear_present_landscapes")
results <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation/results/RDA-2027-14134"
gf_file <- file.path(out, "GF_present_landscape.tif")
rdf_file <- file.path(out, "RDAforest_present_landscape.tif")
env_file <- "/home/nadaline/Documents/Phd/RADseq/Ameiva_Variables_rda.csv"
stopifnot(all(file.exists(c(gf_file, rdf_file, env_file))))

gf <- rast(gf_file)
rdf <- rast(rdf_file)
samples <- read.csv(env_file, check.names = FALSE)
world <- map_data("world")

gf_pca <- readRDS(file.path(out, "GF_present_visualization_PCA.rds"))
rdf_pca <- readRDS(file.path(out, "RDAforest_present_visualization_PCA.rds"))
gf_variance <- 100 * gf_pca$sdev^2 / sum(gf_pca$sdev^2)
rdf_variance <- 100 * rdf_pca$sdev^2 / sum(rdf_pca$sdev^2)

#### 2. DESCRIPTIVE CLIMATIC DIRECTION OF THE VISUALIZATION AXES ####

predictors_gf <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio12", "bio13", "bio15")
predictors_rdf <- c("bio13", "bio3", "bio2")
climate_dir <- "/projects/rstudio-export/Ameiva"
template <- rast(file.path(results, "spatial_offset", "18_RDA_present_axes.tif"))[[1]]
climate_present <- lapply(predictors_gf, function(p) {
  x <- rast(file.path(climate_dir, paste0("CHELSA_", p, "_1981-2010_SA.tif")))
  x <- aggregate(crop(x, ext(-85, -30, -35, 15)), fact = 5, fun = mean, na.rm = TRUE)
  stopifnot(compareGeom(x, template, stopOnError = FALSE))
  mask(x, template)
})
names(climate_present) <- predictors_gf

# Descriptive rank correlations use original climate at exactly the finite
# cells/scores shown by each existing PC1 surface. They are not model tests.
gf_values <- values(gf, mat = FALSE)
gf_direction <- do.call(rbind, lapply(predictors_gf, function(p) {
  climate_values <- values(climate_present[[p]], mat = FALSE)
  keep <- is.finite(gf_values) & is.finite(climate_values)
  rho <- cor(gf_values[keep], climate_values[keep], method = "spearman")
  loading <- gf_pca$rotation[p, 1]
  data.frame(Method = "Gradient Forest", Predictor = toupper(p), rho_with_PC1 = rho, abs_rho = abs(rho),
    PC1_loading = loading, sign_agreement = sign(rho) == sign(loading),
    PC1_end_associated_with_higher_values = ifelse(rho >= 0, "Orange end", "Purple end"))
}))
rdf_df <- as.data.frame(rdf, xy = TRUE, na.rm = TRUE); names(rdf_df)[3] <- "score"
rdf_direction <- do.call(rbind, lapply(predictors_rdf, function(p) {
  climate_values <- extract(climate_present[[p]], rdf_df[, c("x", "y")], ID = FALSE)[, 1]
  keep <- is.finite(rdf_df$score) & is.finite(climate_values)
  rho <- cor(rdf_df$score[keep], climate_values[keep], method = "spearman")
  data.frame(Method = "RDAforest", Predictor = toupper(p), rho_with_PC1 = rho, abs_rho = abs(rho),
    PC1_loading = NA_real_, sign_agreement = NA, PC1_end_associated_with_higher_values = ifelse(rho >= 0, "Orange end", "Purple end"))
}))
gf_direction <- gf_direction[order(-gf_direction$abs_rho), ]
rdf_direction <- rdf_direction[order(-rdf_direction$abs_rho), ]
direction <- rbind(gf_direction, rdf_direction)
direction <- direction[order(direction$Method, -direction$abs_rho), ]
write.csv(direction, file.path(out, "Figure_4_PC1_climate_direction.csv"), row.names = FALSE)

bio_name <- c(BIO1 = "annual mean temperature", BIO2 = "mean diurnal range", BIO3 = "isothermality",
  BIO4 = "temperature seasonality", BIO5 = "maximum temperature of warmest month",
  BIO6 = "minimum temperature of coldest month", BIO7 = "annual temperature range",
  BIO12 = "annual precipitation", BIO13 = "precipitation of wettest month", BIO15 = "precipitation seasonality")
climatic_end <- function(x, end) {
  z <- x[x$PC1_end_associated_with_higher_values == end, , drop = FALSE]
  if (!nrow(z)) return("no retained predictor")
  paste(unname(bio_name[z$Predictor[seq_len(min(3, nrow(z)))]]), collapse = ", ")
}

#### 3. DISPLAY-ONLY CARTOGRAPHY AND CLIMATE-ASSOCIATION INSETS ####

# No displayed value is interpolated or smoothed. GF keeps its native 1/24-degree
# pixels; RDAforest retains its validated native 4,532-cell prediction grid.
map_panel <- function(x, title, legend, display_data = NULL) {
  q <- quantile(values(x, mat = FALSE), c(.02, .98), na.rm = TRUE)
  lim <- max(abs(q)) * c(-1, 1)
  if (is.null(display_data)) {
    d <- as.data.frame(x, xy = TRUE, na.rm = TRUE)
    names(d)[3] <- "score"
  } else d <- display_data
  ggplot(d, aes(x, y, fill = score)) +
    geom_raster() +
    geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE,
      fill = NA, colour = "grey72", linewidth = .16) +
    geom_point(data = samples, aes(long, lat), inherit.aes = FALSE, shape = 21,
      fill = "white", colour = "#5B2A86", stroke = .18, size = .65, alpha = .62) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15), expand = FALSE) +
    scale_fill_gradient2(low = "#32134F", mid = "#FFF7EF", high = "#F28E2B",
      midpoint = 0, limits = lim, oob = scales::squish, name = legend,
      guide = guide_colorbar(barwidth = unit(38, "mm"), barheight = unit(3.5, "mm"))) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", colour = "#32134F", size = 12),
      axis.title = element_text(size = 8), axis.text = element_text(size = 7),
      legend.position = "bottom", legend.title = element_text(face = "bold", size = 8),
      plot.margin = margin(3, 3, 1, 3, "mm"))
}

# RDAforest has 4,532 native prediction cells. Bilinear resampling is used only
# to render the native cells more smoothly; the analytical raster remains the
# surface used above for all correlations and reported statistics.
fine_rdf <- rast(ext(rdf), resolution = .10, crs = crs(rdf))
fine_land <- resample(template, fine_rdf, method = "near")
rdf_display_raster <- mask(resample(rdf, fine_rdf, method = "bilinear"), fine_land)
rdf_display <- as.data.frame(rdf_display_raster, xy = TRUE, na.rm = TRUE); names(rdf_display)[3] <- "score"
writeRaster(rdf_display_raster, file.path(out, "RDAforest_present_landscape_bilinear_display_only.tif"), overwrite = TRUE)

p_gf <- map_panel(gf, "Gradient Forest", "GF landscape PC1") + labs(tag = "a")
p_rdf <- map_panel(rdf, "RDAforest", "RDAforest landscape PC1", rdf_display) + labs(tag = "b")

climate_inset <- function(x) {
  short_names <- c(BIO7 = "Annual temp. range", BIO4 = "Temp. seasonality", BIO2 = "Diurnal temp. range",
    BIO12 = "Annual precipitation", BIO3 = "Isothermality", BIO13 = "Wettest-month precip.")
  x$Predictor <- factor(unname(short_names[x$Predictor]), levels = rev(unname(short_names[x$Predictor])))
  x$direction <- ifelse(x$rho_with_PC1 < 0, "Purple end", "Orange end")
  ggplot(x, aes(Predictor, rho_with_PC1, fill = direction)) + geom_col(width = .62) +
    geom_hline(yintercept = 0, linewidth = .25, colour = "grey45") + coord_flip() +
    scale_fill_manual(values = c("Purple end" = "#5B2A86", "Orange end" = "#F28E2B"), guide = "none") +
    labs(title = "Climatic direction of PC1", x = NULL, y = "Association with landscape PC1 (Spearman rho)") +
    theme_minimal(base_size = 8) + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", colour = "#32134F", size = 9), axis.title.y = element_text(size = 7),
      axis.text = element_text(size = 7), plot.margin = margin(0, 3, 2, 3, "mm"))
}
inset_gf <- climate_inset(gf_direction[seq_len(4), ])
inset_rdf <- climate_inset(rdf_direction)
left <- (p_gf / inset_gf) + plot_layout(heights = c(5, 1.15))
right <- (p_rdf / inset_rdf) + plot_layout(heights = c(5, 1.15))
figure <- (left | right) &
  theme(plot.tag = element_text(face = "bold", colour = "#32134F", size = 11), legend.position = "bottom")

ggsave(file.path(out, "Figure_4_present_nonlinear_landscapes_final.pdf"), figure,
  width = 12, height = 6.4, device = cairo_pdf)
ggsave(file.path(out, "Figure_4_present_nonlinear_landscapes_final.png"), figure,
  width = 12, height = 6.4, dpi = 600, bg = "white")

#### 4. CAPTION, RESULTS TEXT AND RESOLUTION CHECK ####

caption <- paste0(
  "Figure 4. Present-day nonlinear climate-associated genomic landscapes in ",
  "<i>Ameiva ameiva</i>. (a) Primary axis of the present-day genomic-turnover space derived from the cumulative transformations of the controlled Gradient Forest model (53.5% of variation). ",
  "(b) Primary axis of the present-day RDAforest prediction space (32.3% of variation). Insets show descriptive correlations between each displayed PC1 and original present-day climate. Both axes were derived for visualization of the multivariate model outputs, and their signs are arbitrary. Sampled localities are shown as open circles. Gradient Forest PC2 (22.1%) and RDAforest PC2 (18.8%) are shown in the Supplementary Information. Purple and orange indicate opposite ends of each model-specific axis and are not directly comparable in magnitude between models."
)
writeLines(caption, file.path(out, "Figure_4_caption.txt"))

results_text <- c(
  paste0("The GF landscape axis contrasted northern Amazonia and the Guiana Shield, characterized primarily by higher ", climatic_end(gf_direction, "Purple end"), ", with the Southern Andean margin, characterized primarily by higher ", climatic_end(gf_direction, "Orange end"), " (Figure 4a)."),
  paste0("The RDAforest landscape axis contrasted the central-eastern Brazilian interior, associated primarily with higher ", climatic_end(rdf_direction, "Purple end"), ", with northern and western Amazonia and the Andean margin, associated primarily with higher ", climatic_end(rdf_direction, "Orange end"), " (Figure 4b)."),
  "These are descriptive climatic directions of visualization axes, not estimates of fitness, selection, or causal effects."
)
writeLines(results_text, file.path(out, "Figure_4_results_sentences.txt"))

direction_audit <- c("FIGURE 4 PC1 CLIMATIC-DIRECTION AUDIT", "",
  "GF-PC1 loadings were read from GF_present_visualization_PCA.rds. Spearman correlations relate saved GF-PC1 scores to original present-day CHELSA values at common finite cells.",
  "RDAforest-PC1 correlations relate saved PC1 scores to original present-day CHELSA values extracted at validated RDAforest grid-cell centres.",
  "Correlations are descriptive only; no inferential test, model refit, climate recalculation, or sign reversal was performed.",
  "Purple end denotes lower PC1 scores; orange end denotes higher PC1 scores. Axis signs were retained from raw visualization PCA and are arbitrary.", "",
  "Ranked results are in Figure_4_PC1_climate_direction.csv.")
writeLines(direction_audit, file.path(out, "Figure_4_PC1_climate_direction_audit.txt"))

resolution_audit <- c(
  "FIGURE 4 RESOLUTION AUDIT", "",
  "GF source: visualization-only PC1 from the present-day cumulative-transform surface of the 32 controlled whole-panel (14,134-SNP) GF blocks.",
  "RDAforest source: visualization-only PC1 from saved direct present-day predictions of 12_RDAforest_fit_v2610_nonquarter_geoonly_ecdf_v2_pcs8.rds.",
  sprintf("GF native displayed raster: %d finite cells; resolution %.8f x %.8f degrees.", ncell(gf) - global(is.na(gf), "sum", na.rm = TRUE)[1, 1], res(gf)[1], res(gf)[2]),
  sprintf("RDAforest native displayed raster: %d finite cells; resolution %.8f x %.8f degrees.", ncell(rdf) - global(is.na(rdf), "sum", na.rm = TRUE)[1, 1], res(rdf)[1], res(rdf)[2]),
  "No finer native RDAforest present-day prediction raster was found in the current final spatial-output directory; only the validated 4,532-cell prediction grid exists.",
  "GF is rendered at native raster-cell resolution. RDAforest panel b uses display-only bilinear resampling of the 4,532 native prediction cells; the saved analytical RDAforest raster, PC scores, statistics and signs were not modified.",
  "Axis signs were retained from the raw visualization PCA and are arbitrary.",
  sprintf("GF PC1 %.2f%%; GF PC2 %.2f%%; GF PC1+PC2 %.2f%%.", gf_variance[1], gf_variance[2], sum(gf_variance[1:2])),
  sprintf("RDAforest PC1 %.2f%%; RDAforest PC2 %.2f%%; RDAforest PC1+PC2 %.2f%%.", rdf_variance[1], rdf_variance[2], sum(rdf_variance[1:2]))
)
writeLines(resolution_audit, file.path(out, "Figure_4_resolution_audit.txt"))
