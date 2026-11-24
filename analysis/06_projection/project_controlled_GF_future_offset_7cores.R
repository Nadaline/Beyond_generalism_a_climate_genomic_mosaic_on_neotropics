
# 1. Read controlled Gradient Forest curves, climate data and raster template ----

library(terra)
library(parallel)

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"
results <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation/results/RDA-2027-14134"
source_out <- file.path(project, "results", "nonlinear_present_landscapes")
out <- file.path(results, "spatial_offset", "GF_controlled_wholepanel_future_offset")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

curve_file <- file.path(source_out, "GF_controlled_mean_cumulative_curves.csv")
template_file <- file.path(results, "spatial_offset", "18_RDA_present_axes.tif")
environment_file <- "/home/nadaline/Documents/Phd/RADseq/Ameiva_Variables_rda.csv"
climate_dir <- "/projects/rstudio-export/Ameiva"
stopifnot(all(file.exists(c(curve_file, template_file, environment_file))))

predictors <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio12", "bio13", "bio15")
periods <- c("1981-2010", "2041-2070", "2071-2100")
curves <- read.csv(curve_file, stringsAsFactors = FALSE)
environment <- read.csv(environment_file, check.names = FALSE, stringsAsFactors = FALSE)
template <- rast(template_file)[[1]]
workers <- min(7L, detectCores(logical = TRUE), length(predictors))

lineage <- c(
  "CONTROLLED WHOLE-PANEL GF FUTURE OFFSET", "",
  "Model source: 32 completed controlled 14,134-SNP GF block objects under gradient_forest_all_14134/controlled_PC1_PC3_lon_lat/.",
  "Transformation source: mean cumulative-importance curves saved in results/nonlinear_present_landscapes/GF_controlled_mean_cumulative_curves.csv.",
  "Predictors: BIO1–BIO7, BIO12, BIO13 and BIO15; ECDFs estimated from the original 245 sampled environmental records.",
  "Future climate source: existing CHELSA 2041–2070 and 2071–2100 rasters in /projects/rstudio-export/Ameiva.",
  paste0("Workers requested/used: up to 7 / ", workers, "."),
  "No GF fit, candidate selection, climate generation, RDAforest procedure, or cross-validation is run by this script."
)
writeLines(lineage, file.path(out, "GF_controlled_future_offset_lineage.txt"))

# 2. Transform present and future climatic rasters ----

# Each worker transforms one predictor across the three climatic periods.
# Aggregate the CHELSA layers by five and align them to the reference grid.
transform_one_predictor <- function(predictor) {
  observed <- environment[[paste0("X.CHELSA_", predictor, "_1981.2010_SA")]]
  observed <- sort(observed[is.finite(observed)])
  probability <- seq_along(observed) / length(observed)
  curve <- curves[curves$predictor == predictor, ]
  span <- diff(range(observed))
  output <- character(length(periods))

  for (i in seq_along(periods)) {
    period <- periods[i]
    raw <- rast(file.path(climate_dir, paste0("CHELSA_", predictor, "_", period, "_SA.tif")))
    raw <- crop(raw, ext(-85, -30, -35, 15))
    raw <- aggregate(raw, fact = 5, fun = mean, na.rm = TRUE)
    stopifnot(compareGeom(raw, template, stopOnError = FALSE))
    raw <- mask(raw, template)

    # Allow 10% of the observed predictor range beyond present-day limits and 20%
    # for future climate. Clamp values to curve endpoints and exclude cells beyond
    # the allowed range from the offset calculation.
    buffer <- if (period == "1981-2010") .10 else .20
    outside_allowed <- raw < (min(observed) - buffer * span) | raw > (max(observed) + buffer * span)
    ecdf_raster <- app(raw, function(v) approx(observed, probability, xout = v, rule = 2, ties = "ordered")$y)
    transformed <- app(ecdf_raster, function(v) approx(curve$ecdf, curve$mean_cumimp, xout = v, rule = 2)$y)
    names(transformed) <- predictor
    transformed_file <- file.path(out, paste0("GF_controlled_", predictor, "_", period, "_transformed.tif"))
    status_file <- file.path(out, paste0("GF_controlled_", predictor, "_", period, "_outside_allowed.tif"))
    writeRaster(transformed, transformed_file, overwrite = TRUE)
    writeRaster(outside_allowed, status_file, overwrite = TRUE)
    output[i] <- transformed_file
  }
  output
}

transformed_files <- mclapply(predictors, transform_one_predictor, mc.cores = workers, mc.preschedule = FALSE)
names(transformed_files) <- predictors
stopifnot(all(vapply(transformed_files, function(x) all(file.exists(x)), logical(1))))

# 3. Calculate genomic offset within the supported climatic domain ----

stack_for_period <- function(period) {
  rast(vapply(transformed_files, function(x) x[which(periods == period)], character(1)))
}
present <- stack_for_period("1981-2010")
mid <- stack_for_period("2041-2070")
late <- stack_for_period("2071-2100")
names(present) <- names(mid) <- names(late) <- predictors

allowed_for_period <- function(period) {
  status_files <- file.path(out, paste0("GF_controlled_", predictors, "_", period, "_outside_allowed.tif"))
  ifel(sum(rast(status_files)) == 0, 1, NA)
}
allowed_present <- allowed_for_period("1981-2010")
allowed_mid <- allowed_for_period("2041-2070")
allowed_late <- allowed_for_period("2071-2100")
allowed_mid_combined <- allowed_present * allowed_mid
allowed_late_combined <- allowed_present * allowed_late

# Euclidean distance in cumulative-importance space, masked in both periods.
offset_mid <- sqrt(sum((mid - present)^2)) * allowed_mid_combined
offset_late <- sqrt(sum((late - present)^2)) * allowed_late_combined
names(offset_mid) <- "GF_controlled_wholepanel_offset_2041_2070"
names(offset_late) <- "GF_controlled_wholepanel_offset_2071_2100"
writeRaster(offset_mid, file.path(out, "GF_controlled_wholepanel_offset_2041_2070.tif"), overwrite = TRUE)
writeRaster(offset_late, file.path(out, "GF_controlled_wholepanel_offset_2071_2100.tif"), overwrite = TRUE)
writeRaster(allowed_mid_combined, file.path(out, "GF_controlled_projection_domain_2041_2070.tif"), overwrite = TRUE)
writeRaster(allowed_late_combined, file.path(out, "GF_controlled_projection_domain_2071_2100.tif"), overwrite = TRUE)

# 4. Summarise domain coverage and offset distributions ----
domain_summary <- do.call(rbind, Map(function(period, x) {
  values_x <- values(x, mat = FALSE)
  data.frame(period = period, terrestrial_cells = length(values_x), supported_cells = sum(is.finite(values_x)), excluded_cells = sum(!is.finite(values_x)), excluded_proportion = mean(!is.finite(values_x)))
}, c("2041-2070", "2071-2100"), list(allowed_mid_combined, allowed_late_combined)))
write.csv(domain_summary, file.path(out, "GF_controlled_projection_domain_summary.csv"), row.names = FALSE)

offset_summary <- do.call(rbind, Map(function(period, x) {
  v <- values(x, mat = FALSE)
  data.frame(period = period, finite_cells = sum(is.finite(v)), q10 = quantile(v, .10, na.rm = TRUE), median = median(v, na.rm = TRUE), q90 = quantile(v, .90, na.rm = TRUE))
}, c("2041-2070", "2071-2100"), list(offset_mid, offset_late)))
write.csv(offset_summary, file.path(out, "GF_controlled_wholepanel_offset_summary.csv"), row.names = FALSE)

cat("Completed controlled whole-panel GF future offset post-processing\n")
print(domain_summary)
print(offset_summary)
