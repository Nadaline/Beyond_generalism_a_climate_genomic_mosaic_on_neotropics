# Landscape genomics of Ameiva ameiva
# Main association panel: 245 individuals and 14,134 SNPs.


# 0. Load packages, define paths and set analysis options ----

.libPaths(c(normalizePath("R-library"), .libPaths()))

library(vcfR)
library(dartR.base)
library(dartR.popgen)
library(dartR.spatial)
library(adegenet)
library(vegan)
library(LEA)
library(Boruta)
library(usdm)
library(qvalue)
library(robust)
library(ggplot2)
library(cluster)
library(terra)
library(gradientForest)
library(RDAforest)
library(GenomicRanges)
library(IRanges)
library(maps)
library(sf)

ameiva_world_map <- ggplot2::map_data("world")

set.seed(2027)

idv_purple_dark <- "#2F1647"
idv_purple <- "#5B2A86"
idv_purple_light <- "#B99AD0"
idv_orange_dark <- "#C76100"
idv_orange <- "#F28E2B"
idv_orange_light <- "#F8C58F"
idv_teal <- "#126E5F"
idv_cream <- "#FFF7EF"
idv_grey <- "#BDBDBD"

idv_sequential <- c(
  idv_cream, idv_orange_light, idv_orange,
  "#A95572", idv_purple, idv_purple_dark
)
idv_diverging <- c(
  idv_purple_dark, idv_purple, idv_purple_light,
  "#FAF7F2",
  idv_orange_light, idv_orange, idv_orange_dark
)

idv_discrete <- function(n) {
  anchors <- c(
    idv_purple_dark, idv_purple, idv_purple_light,
    idv_orange_light, idv_orange, idv_orange_dark
  )
  if (n == 1) return(idv_purple)
  if (n == 2) return(c(idv_purple, idv_orange))
  if (n == 3) return(c(idv_purple, idv_orange, idv_teal))
  if (n == 4) return(c(
    idv_purple_dark, idv_purple_light,
    idv_orange_light, idv_orange_dark
  ))
  if (n <= length(anchors)) return(anchors[seq_len(n)])
  grDevices::colorRampPalette(anchors)(n)
}

theme_ameiva <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey91", linewidth = .3),
      plot.title = element_text(face = "bold", colour = idv_purple_dark),
      plot.subtitle = element_text(colour = "grey35"),
      strip.text = element_text(face = "bold", colour = idv_purple_dark),
      strip.background = element_rect(
        fill = "#F3EAF8", colour = NA
      ),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(colour = "grey20")
    )
}

raster_to_ggplot <- function(x, value_name = "value", period = NULL) {
  result <- as.data.frame(x, xy = TRUE, na.rm = FALSE)
  names(result)[3] <- value_name
  if (!is.null(period)) result$period <- period
  result
}

# Encode three ordination axes in a shared false-colour scale.
idv_false_colour <- function(red, green, blue) {
  channels <- cbind(red, green, blue)
  channels[!is.finite(channels)] <- 0
  channels[channels < 0] <- 0
  channels[channels > 1] <- 1
  anchors <- t(grDevices::col2rgb(c(
    "#6A3D9A", # purple: axis 1
    "#E69F00", # orange: axis 2 (Okabe-Ito)
    "#00796B"  # dark teal: axis 3
  )) / 255)
  weights <- channels^1.7 + .015
  mixed <- weights %*% anchors / rowSums(weights)
  grDevices::rgb(mixed[, 1], mixed[, 2], mixed[, 3])
}

adaptation_ggplot_panels <- function(adaptation, environment) {
  pc <- as.data.frame(adaptation$PCs[, 1:2, drop = FALSE])
  names(pc) <- c("PC1", "PC2")
  xy <- as.data.frame(adaptation$xy)
  names(xy)[1:2] <- c("longitude", "latitude")
  clustered <- length(adaptation$clusters) > 1

  fitted_environment <- vegan::envfit(
    pc, environment, permutations = 0
  )
  arrows <- as.data.frame(vegan::scores(
    fitted_environment, display = "vectors", choices = 1:2
  ))
  names(arrows)[1:2] <- c("PC1", "PC2")
  arrows$variable <- rownames(arrows)
  arrow_multiplier <- .60 * min(
    max(abs(pc$PC1), na.rm = TRUE) /
      max(abs(arrows$PC1), na.rm = TRUE),
    max(abs(pc$PC2), na.rm = TRUE) /
      max(abs(arrows$PC2), na.rm = TRUE)
  )
  arrows$PC1 <- arrows$PC1 * arrow_multiplier
  arrows$PC2 <- arrows$PC2 * arrow_multiplier

  if (clustered) {
    pc$group <- factor(adaptation$clusters)
    xy$group <- pc$group
    group_colours <- setNames(
      idv_discrete(nlevels(pc$group)), levels(pc$group)
    )
    ordination_plot <- ggplot(pc, aes(PC1, PC2, colour = group)) +
      geom_point(size = 1.7, alpha = .78) +
      scale_colour_manual(values = group_colours, guide = "none")
    map_plot <- ggplot(xy, aes(longitude, latitude, colour = group)) +
      geom_polygon(data = ameiva_world_map, aes(long, lat, group = group), inherit.aes = FALSE, fill = idv_cream, colour = "grey78", linewidth = .2) +
      geom_point(shape = 15, size = 1.25, alpha = .88) +
      scale_colour_manual(values = group_colours, guide = "none")
    cluster_sizes <- as.data.frame(table(group = pc$group))
    size_plot <- ggplot(
      cluster_sizes, aes(group, Freq, fill = group)
    ) +
      geom_col(width = .72) +
      scale_fill_manual(values = group_colours, guide = "none") +
      theme_ameiva() +
      labs(title = "Neighborhood sizes", x = "Cluster", y = "Cells")
  } else {
    # Retain information from all three axes in the colour mapping.
    original_channels <- t(grDevices::col2rgb(adaptation$colors)) / 255
    pc$brand_colour <- idv_false_colour(
      original_channels[, 1],
      original_channels[, 2],
      original_channels[, 3]
    )
    xy$brand_colour <- pc$brand_colour
    ordination_plot <- ggplot(
      pc, aes(PC1, PC2, colour = brand_colour)
    ) +
      geom_point(size = 1.7, alpha = .78) +
      scale_colour_identity()
    grid_dx <- median(diff(sort(unique(xy$longitude))))
    grid_dy <- median(diff(sort(unique(xy$latitude))))
    map_plot <- ggplot(
      xy, aes(longitude, latitude, fill = brand_colour)
    ) +
      geom_polygon(data = ameiva_world_map, aes(long, lat, group = group), inherit.aes = FALSE, fill = idv_cream, colour = "grey78", linewidth = .2) +
      geom_tile(
        width = grid_dx * 1.015, height = grid_dy * 1.015,
        linewidth = 0, show.legend = FALSE
      ) +
      scale_fill_identity()
    size_plot <- NULL
  }

  ordination_plot <- ordination_plot +
    geom_segment(
      data = arrows,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(2.2, "mm")),
      colour = idv_orange, linewidth = .75
    ) +
    geom_text(
      data = arrows,
      aes(PC1, PC2, label = variable),
      inherit.aes = FALSE,
      colour = idv_orange_dark, fontface = "bold",
      size = 3, nudge_y = .025 * diff(range(pc$PC2))
    ) +
    coord_equal() +
    theme_ameiva() +
    labs(title = "Genomic space", colour = "Cluster")

  map_plot <- map_plot +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Geographic space",
      x = "Longitude", y = "Latitude", colour = "Cluster"
    )

  Filter(Negate(is.null), list(
    ordination = ordination_plot,
    map = map_plot,
    sizes = size_plot
  ))
}

save_ggplot_row <- function(plots, filename, width = 2400,
                            height = 1200, res = 220) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  number_plots <- length(plots)
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = grid::viewport(
      x = (i - .5) / number_plots,
      y = .5,
      width = .98 / number_plots,
      height = .96
    ))
  }
  dev.off()
}

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"
ameiva <- "/home/nadaline/Documents/Phd/RADseq/Ameiva"
hd <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation"

raw_vcf <- file.path(ameiva, "Ameiva_genome", "populations.snps.vcf")
environment_file <- "/home/nadaline/Documents/Phd/RADseq/Ameiva_Variables_rda.csv"
metadata_file <- file.path(ameiva, "Ameiva_Metadata.csv")

out <- file.path(hd, "results", "RDA-2027-14134")
figures <- file.path(project, "results", "figures", "RDA-2027-14134")

dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

# Save the output of long-running analyses in a log file.
execution_log <- file.path(out, "00_execution.log")
execution_connection <- file(execution_log, open = "at")
sink(execution_connection, split = FALSE)
sink(execution_connection, type = "message")

stopifnot(file.exists(raw_vcf))
stopifnot(file.exists(environment_file))
stopifnot(file.exists(metadata_file))
stopifnot(dir.exists(hd))

# Saved intermediate objects allow the analysis to resume. Environmental
# model files are identified by the predictor-selection specification.
USE_SAVED_RDS <- TRUE
# The partial RDA locus scan uses a family-wise Bonferroni threshold of 0.01.
ENVIRONMENT_VERSION <- "nonquarter_geoonly_bonferroni_v1"

# Select which analytical modules to run.
RUN_SNMF <- TRUE
PREPARE_TREEMIX <- TRUE
PREPARE_STAIRWAYPLOT2 <- TRUE
RUN_LFMM <- TRUE
RUN_GRADIENT_FOREST <- TRUE
RUN_RDAFOREST <- TRUE
RUN_RDAFOREST_SPATIAL <- TRUE
RUN_PARTITIONED_RDA_TESTS <- FALSE
RUN_SPATIAL_PROJECTION <- TRUE
SNMF_MINIMUM_REPETITIONS <- 10
SNMF_SELECTED_K <- 4
SNMF_K_VALUES <- 1:8
RDAFOREST_PCS_OVERRIDE <- suppressWarnings(as.integer(
  Sys.getenv("RDAFOREST_PCS_OVERRIDE", unset = "0")
))


# 1. Read environmental data and individual metadata ----

Variables <- read.csv(
  environment_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(Variables) <- Variables$InternalID

Ameiva_metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(all(c("id", "pop") %in% names(Ameiva_metadata)))
stopifnot(!anyDuplicated(Ameiva_metadata$id))



# 2. Filter the VCF and match individuals to climate data ----

# Apply locus and individual call-rate filters before removing invariant loci.
dart_cache_file <- file.path(out, "02_genlight_after_dartR_filters.rds")

if (USE_SAVED_RDS && file.exists(dart_cache_file)) {
  ref_dart <- readRDS(dart_cache_file)
} else {
  ref_vcf <- read.vcfR(raw_vcf, verbose = FALSE)
  ref_dart <- vcfR2genlight(ref_vcf)
  ploidy(ref_dart) <- 2
  ref_dart <- gl.compliance.check(ref_dart)
  ref_dart <- gl.add.indmetrics(ref_dart, metadata_file, verbose = 0)

  ref_dart <- gl.filter.callrate(
    ref_dart, method = "loc", threshold = .70, verbose = 0
  )
  ref_dart <- gl.filter.callrate(
    ref_dart, method = "ind", threshold = .70, verbose = 0
  )
  ref_dart <- gl.filter.monomorphs(ref_dart, verbose = 0)
  ref_dart <- gl.filter.allna(ref_dart, verbose = 0)

  saveRDS(ref_dart, dart_cache_file)
}

nInd(ref_dart)                      # 250 before matching climate data
nLoc(ref_dart)

# Five genetically filtered individuals lack corresponding environmental data.
ids <- indNames(ref_dart)
ids <- ids[ids %in% Variables$InternalID]

length(ids)                         # 245 matched individuals
stopifnot(length(ids) == 245)
stopifnot(!anyDuplicated(ids))

ref_dart <- gl.keep.ind(ref_dart, ind.list = ids, verbose = 0)
Variables <- Variables[ids, , drop = FALSE]
stopifnot(identical(indNames(ref_dart), rownames(Variables)))

sample_file <- file.path(out, "01_individuals_245.txt")
writeLines(ids, sample_file)

write.csv(
  data.frame(
    InternalID = ids,
    population = Variables$POP
  ),
  file.path(out, "01_individual_registry.csv"),
  row.names = FALSE
)


# 3.1 Filter SNPs within the matched set of individuals ----

# Recalculate missingness and allele counts after subsetting individuals.
# Retain loci with at least three copies of each allele and at most 30% missing data.

filtered_vcf <- file.path(out, "02_Ameiva_245_missing30_MAC3.vcf.gz")

if (!USE_SAVED_RDS || !file.exists(filtered_vcf)) {

  filter_command <- paste(
    "bcftools view",
    "-S", shQuote(sample_file),
    "-Ou", shQuote(raw_vcf),
    "| bcftools +fill-tags -Ou -- -t AC,AN,F_MISSING",
    "| bcftools view -i",
    shQuote("F_MISSING<=0.30 && AC>=3 && (AN-AC)>=3"),
    "-Oz -o", shQuote(filtered_vcf)
  )

  status_filter <- system(filter_command)
  stopifnot(status_filter == 0)

  # The Stacks VCF contains unsorted scaffolds; vcfR reads it sequentially.
}

number_snps <- as.integer(system2(
  "bcftools",
  c("view", "-H", shQuote(filtered_vcf)),
  stdout = TRUE
) |> length())

number_snps                         # 14,134 SNPs
stopifnot(number_snps == 14134)


# 3.2 Convert VCF genotypes to an individual-by-SNP matrix ----

genotype_cache_file <- file.path(out, "03_genotype_matrix_14134.rds")
reuse_genotype_cache <- USE_SAVED_RDS && file.exists(genotype_cache_file)

if (reuse_genotype_cache) {
  genotype_cache <- readRDS(genotype_cache_file)
  geno <- genotype_cache$geno
  variant_table <- genotype_cache$variant_table
  reuse_genotype_cache <-
    identical(rownames(geno), ids) &&
    nrow(geno) == 245 && ncol(geno) == number_snps &&
    nrow(variant_table) == number_snps
}

if (!reuse_genotype_cache) {

  vcf <- read.vcfR(filtered_vcf, verbose = FALSE)

  genotype_strings <- extract.gt(
    vcf,
    element = "GT",
    as.numeric = FALSE,
    return.alleles = FALSE
  )

  # Code diploid genotypes as counts of the alternative allele (0, 1 or 2).
  geno <- matrix(
    NA_real_,
    nrow = nrow(genotype_strings),
    ncol = ncol(genotype_strings),
    dimnames = dimnames(genotype_strings)
  )
  geno[genotype_strings %in% c("0/0", "0|0")] <- 0
  geno[genotype_strings %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  geno[genotype_strings %in% c("1/1", "1|1")] <- 2

  unexpected_GT <- setdiff(
    unique(as.vector(genotype_strings)),
    c("0/0", "0|0", "0/1", "1/0", "0|1", "1|0", "1/1", "1|1", NA, "./.", ".|.")
  )
  stopifnot(length(unexpected_GT) == 0)

  geno <- t(geno)
  storage.mode(geno) <- "numeric"

  variant_table <- as.data.frame(getFIX(vcf), stringsAsFactors = FALSE)
  colnames(geno) <- variant_table$ID

  stopifnot(all(ids %in% rownames(geno)))
  geno <- geno[ids, , drop = FALSE]

  saveRDS(
    list(geno = geno, variant_table = variant_table),
    genotype_cache_file
  )

  rm(vcf, genotype_strings)
  gc()
}

stopifnot(nrow(geno) == 245)
stopifnot(ncol(geno) == 14134)
stopifnot(identical(rownames(geno), rownames(Variables)))


# 4. Summarise missingness, allele frequencies and RAD loci ----

missing_locus <- colMeans(is.na(geno))
called_alleles <- 2 * colSums(!is.na(geno))
alt_count <- colSums(geno, na.rm = TRUE)
minor_allele_count <- pmin(alt_count, called_alleles - alt_count)
allele_frequency <- alt_count / called_alleles
maf <- pmin(allele_frequency, 1 - allele_frequency)
rad_locus <- sub(":.*", "", colnames(geno))

qc_loci <- data.frame(
  stacks_id = colnames(geno),
  scaffold = variant_table$CHROM,
  position = as.integer(variant_table$POS),
  rad_locus = rad_locus,
  missing = missing_locus,
  MAC = minor_allele_count,
  MAF = maf
)

write.csv(qc_loci, file.path(out, "02_panel_14134_QC.csv"), row.names = FALSE)

length(unique(rad_locus))           # approximately 5,776 RAD loci
summary(missing_locus)
summary(minor_allele_count)
summary(rowMeans(is.na(geno)))

stopifnot(max(missing_locus) <= .300001)
stopifnot(min(minor_allele_count) >= 3)


# 5. Select one SNP per RAD locus for linkage sensitivity ----

# Compare genome-wide PCA with a panel that retains one SNP per RAD tag.
# Within each tag, favour lower missingness and then higher minor allele count.

one_per_RAD <- qc_loci[order(
  qc_loci$rad_locus,
  qc_loci$missing,
  -qc_loci$MAC
),]

one_per_RAD <- one_per_RAD[!duplicated(one_per_RAD$rad_locus),]
geno_unlinked <- geno[, one_per_RAD$stacks_id, drop = FALSE]

ncol(geno_unlinked)                 # approximately 5,776 SNPs

write.csv(
  one_per_RAD,
  file.path(out, "02_panel_one_SNP_per_RAD_QC.csv"),
  row.names = FALSE
)


# 6. Estimate population structure and impute missing genotypes with sNMF ----

# Impute missing genotypes using the fitted sNMF ancestry model.

snmf_file <- file.path(out, "03_Ameiva_14134.lfmm")
snmf_matrix <- geno
snmf_matrix[is.na(snmf_matrix)] <- 9
LEA::write.lfmm(snmf_matrix, snmf_file)

snmf_project_file <- sub("\\.lfmm$", ".snmfProject", snmf_file)

if (RUN_SNMF) {

  set.seed(2027)

  if (USE_SAVED_RDS && file.exists(snmf_project_file)) {
    snmf_project <- load.snmfProject(snmf_project_file)
    existing_runs <- min(sapply(SNMF_K_VALUES, function(k) {
      length(cross.entropy(snmf_project, K = k))
    }))
  } else {
    existing_runs <- 0
  }

  if (existing_runs < SNMF_MINIMUM_REPETITIONS) {
    snmf_project <- snmf(
      snmf_file,
      K = SNMF_K_VALUES,
      repetitions = SNMF_MINIMUM_REPETITIONS - existing_runs,
      entropy = TRUE,
      project = if (existing_runs > 0) "continue" else "new",
      CPU = 6,
      seed = 2027
    )
  }

  cross_entropy <- do.call(rbind, lapply(SNMF_K_VALUES, function(k) {
    entropy_k <- as.numeric(cross.entropy(snmf_project, K = k))
    data.frame(
      K = k,
      run = seq_along(entropy_k),
      cross_entropy = entropy_k
    )
  }))

  cross_entropy_summary <- aggregate(
    cross_entropy ~ K,
    cross_entropy,
    function(x) c(mean = mean(x), minimum = min(x), sd = sd(x))
  )

  cross_entropy_summary <- data.frame(
    K = cross_entropy_summary$K,
    mean = cross_entropy_summary$cross_entropy[, "mean"],
    minimum = cross_entropy_summary$cross_entropy[, "minimum"],
    sd = cross_entropy_summary$cross_entropy[, "sd"]
  )

  # Mean cross-entropy reaches its minimum at the upper end of the tested K.
  # K = 4 describes the broad structure; both choices are recorded for comparison.
  numerical_minimum_K <- cross_entropy_summary$K[
    cross_entropy_summary$K >= 2
  ][which.min(cross_entropy_summary$mean[cross_entropy_summary$K >= 2])]
  bestK <- SNMF_SELECTED_K
  stopifnot(bestK %in% cross_entropy_summary$K)

  write.csv(
    data.frame(
      selected_K = bestK,
      numerical_minimum_mean_cross_entropy_K = numerical_minimum_K,
      rationale = "K=4: inflection, interpretability and previous validated analysis"
    ),
    file.path(out, "03_sNMF_K_decision.csv"),
    row.names = FALSE
  )

  best_run <- which.min(cross.entropy(snmf_project, K = bestK))

  Q <- LEA::Q(snmf_project, K = bestK, run = best_run)
  # Normalise ancestry coefficients to sum to one for each individual.
  Q <- sweep(Q, 1, rowSums(Q), "/")
  rownames(Q) <- rownames(geno)

  write.csv(cross_entropy, file.path(out, "03_sNMF_cross_entropy_runs.csv"), row.names = FALSE)
  write.csv(cross_entropy_summary, file.path(out, "03_sNMF_cross_entropy_summary.csv"), row.names = FALSE)
  write.csv(cbind(InternalID = rownames(Q), Q), file.path(out, "03_sNMF_Q.csv"), row.names = FALSE)

  plot_entropy <- ggplot(cross_entropy_summary, aes(K, mean)) +
    geom_ribbon(
      aes(ymin = mean - sd, ymax = mean + sd),
      fill = idv_orange_light, alpha = .55
    ) +
    geom_line(linewidth = 1, colour = idv_purple) +
    geom_point(size = 2.6, colour = idv_purple_dark) +
    geom_vline(
      xintercept = bestK, linetype = 2,
      colour = idv_orange, linewidth = .9
    ) +
    scale_x_continuous(breaks = 1:10) +
    theme_ameiva() +
    labs(
      title = "sNMF model selection",
      subtitle = paste("Biologically selected K =", bestK, "| numerical minimum =", numerical_minimum_K),
      y = "Mean cross-entropy"
    )

  ggsave(file.path(figures, "01_sNMF_cross_entropy.png"), plot_entropy,
         width = 7, height = 5, dpi = 400, bg = "white")

  # Use modal genotypes conditional on the fitted ancestry model.
  impute(
    snmf_project,
    snmf_file,
    method = "mode",
    K = bestK,
    run = best_run
  )

  imputed_file <- paste0(snmf_file, "_imputed.lfmm")
  geno_complete <- as.matrix(read.lfmm(imputed_file))
  rownames(geno_complete) <- rownames(geno)
  colnames(geno_complete) <- colnames(geno)

  stopifnot(!anyNA(geno_complete))
}


# 7. Compare genome-wide PCA and the one-SNP-per-RAD panel ----

# PCA is an unconstrained RDA on the scaled, imputed allele-dosage matrix.

pca_genome <- rda(geno_complete, scale = TRUE)

PCA_scores <- as.data.frame(scores(
  pca_genome,
  display = "sites",
  choices = 1:10,
  scaling = 0
))

PCA_scores$InternalID <- rownames(PCA_scores)
PCA_scores$population <- Variables[PCA_scores$InternalID, "POP"]
PCA_scores$ancestry_cluster <- factor(max.col(Q))

write.csv(PCA_scores, file.path(out, "04_PCA_14134_scores.csv"), row.names = FALSE)
saveRDS(pca_genome, file.path(out, "04_PCA_14134.rds"))

plot_pca <- ggplot(PCA_scores, aes(PC1, PC2, colour = ancestry_cluster)) +
  geom_point(size = 2.6, alpha = .9) +
  stat_ellipse(aes(group = ancestry_cluster), linewidth = .5, level = .80) +
  scale_colour_manual(values = idv_discrete(bestK)) +
  theme_ameiva() +
  labs(
    title = "Genomic structure across Ameiva ameiva",
    subtitle = "14,134 SNPs; colours are dominant sNMF ancestry components",
    colour = "sNMF cluster"
  )

ggsave(file.path(figures, "02_PCA_sNMF_structure.png"), plot_pca,
       width = 8, height = 6, dpi = 400, bg = "white")

# Map the predominant genetic cluster and retain the continuous ancestry
# coefficients in the individual bar plot.
structure_map_data <- data.frame(
  longitude = Variables$long,
  latitude = Variables$lat,
  ancestry_cluster = factor(max.col(Q))
)
plot_population_map <- ggplot(
  structure_map_data,
  aes(longitude, latitude, fill = ancestry_cluster)
) +
  geom_polygon(data = ameiva_world_map, aes(long, lat, group = group), inherit.aes = FALSE, fill = idv_cream, colour = "grey75", linewidth = .2) +
  geom_point(
    shape = 21, colour = "white", stroke = .35, size = 2.7,
    position = position_jitter(width = .10, height = .10, seed = 2027)
  ) +
  scale_fill_manual(values = idv_discrete(bestK)) +
  coord_quickmap(xlim = c(-82, -32), ylim = c(-35, 12)) +
  theme_ameiva() +
  theme(panel.grid = element_blank()) +
  labs(
    title = "Continental sampling and ancestry",
    x = "Longitude", y = "Latitude", fill = "Dominant sNMF component"
  )

structure_order <- order(max.col(Q), Variables$long)
Q_long <- do.call(rbind, lapply(seq_len(bestK), function(component) {
  data.frame(
    sample_order = seq_along(structure_order),
    ancestry_component = factor(component, levels = seq_len(bestK)),
    ancestry = Q[structure_order, component]
  )
}))
plot_snmf_ancestry <- ggplot(
  Q_long,
  aes(sample_order, ancestry, fill = ancestry_component)
) +
  geom_col(width = 1.01, colour = NA) +
  scale_fill_manual(values = idv_discrete(bestK)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_ameiva() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("Continuous ancestry at sNMF K = ", bestK),
    x = "Individuals ordered by dominant component and longitude",
    y = "Ancestry proportion", fill = "Component"
  )

ggsave(
  file.path(figures, "32_population_sampling_map.png"),
  plot_population_map, width = 7, height = 6,
  dpi = 400, bg = "white"
)
ggsave(
  file.path(figures, "33_sNMF_ancestry_barplot.png"),
  plot_snmf_ancestry, width = 11, height = 4.2,
  dpi = 400, bg = "white"
)

# Repeat PCA using one SNP per RAD locus.
geno_unlinked_complete <- geno_complete[, colnames(geno_unlinked), drop = FALSE]
pca_unlinked <- rda(geno_unlinked_complete, scale = TRUE)

PCA_unlinked_scores <- as.data.frame(scores(
  pca_unlinked,
  display = "sites",
  choices = 1:3,
  scaling = 0
))

PCA_unlinked_scores$InternalID <- rownames(PCA_unlinked_scores)

write.csv(PCA_unlinked_scores, file.path(out, "04_PCA_one_per_RAD_scores.csv"), row.names = FALSE)
saveRDS(pca_unlinked, file.path(out, "04_PCA_one_per_RAD.rds"))

PCA_concordance <- data.frame(
  axis = paste0("PC", 1:3),
  correlation = sapply(1:3, function(j) {
    abs(cor(PCA_scores[[paste0("PC", j)]], PCA_unlinked_scores[[paste0("PC", j)]]))
  })
)

write.csv(PCA_concordance, file.path(out, "04_PCA_linkage_sensitivity.csv"), row.names = FALSE)


# 8. Prepare geographic group allele counts for TreeMix ----

# TreeMix uses geographic groups from the original metadata rather than sNMF
# assignments. The Parapr and Cerado labels are combined into Cerrado.
# The exported allele counts are the input for the separate TreeMix workflow.
population_treemix <- Ameiva_metadata$pop[
  match(rownames(geno), Ameiva_metadata$id)
]
population_treemix[population_treemix %in% c("Parapr", "Cerado")] <-
  "Cerrado"
population <- factor(population_treemix)
stopifnot(!anyNA(population))

treemix_population_map <- data.frame(
  InternalID = rownames(geno),
  original_population = Ameiva_metadata$pop[
    match(rownames(geno), Ameiva_metadata$id)
  ],
  TreeMix_population = as.character(population)
)
write.csv(
  treemix_population_map,
  file.path(out, "05_TreeMix_population_map.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(TreeMix_population = population)),
  file.path(out, "05_TreeMix_population_sizes.csv"),
  row.names = FALSE
)

if (PREPARE_TREEMIX) {

  # Retain missing genotypes so allele counts use observed copies only.
  gl_treemix <- adegenet::as.genlight(geno)
  adegenet::pop(gl_treemix) <- population
  gl_treemix@ind.names <- rownames(geno)
  gl_treemix@loc.names <- colnames(geno)

  # Initialise the dartR metrics required to calculate allele frequencies.
  gl_treemix <- dartR::gl.compliance.check(gl_treemix, verbose = 0)
  gl_treemix <- dartR::gl.recalc.metrics(gl_treemix, verbose = 0)

  dartR::gl2treemix(
    gl_treemix,
    outfile = "05_Ameiva_14134_treemix.gz",
    outpath = out,
    verbose = 2
  )

  treemix_file <- file.path(out, "05_Ameiva_14134_treemix.gz")
  stopifnot(file.exists(treemix_file), file.info(treemix_file)$size > 0)

  writeLines(
    c(
      "Entrada preparada por dartR::gl2treemix().",
      paste("Individuos:", nrow(geno)),
      paste("SNPs:", ncol(geno)),
      paste("Populacoes:", nlevels(population)),
      "Parapr e Cerado foram reunidas como Cerrado.",
      paste("Maximo recomendado de migracoes:", nlevels(population) - 2),
      paste("Arquivo:", treemix_file),
      paste("Bash disponivel:", file.path(ameiva, "run_treemix_pipeline.sh")),
      "O Bash nao e executado pelo RDA-2027.R."
    ),
    file.path(out, "05_TreeMix_README.txt")
  )
}


# 9. Construct folded site-frequency spectra for Stairway Plot 2 ----

# Use genotypes before the MAC filter to retain rare variants informative
# about recent demographic history. Group individuals by dominant sNMF ancestry.
# Absolute scales depend on mutation rate, generation time and assumed callable RAD length.
STAIRWAY_MUTATION_RATE <- 1.0e-8
STAIRWAY_GENERATION_YEARS <- 2
STAIRWAY_RAD_TAG_LENGTH <- 45
STAIRWAY_TRAINING_PROPORTION <- .67
STAIRWAY_BOOTSTRAPS <- 200

if (PREPARE_STAIRWAYPLOT2) {

  stairway_out <- file.path(out, "stairwayplot2")
  dir.create(stairway_out, recursive = TRUE, showWarnings = FALSE)

  stairway_geno <- as.matrix(ref_dart)
  stairway_population <- paste0(
    "K", max.col(Q[rownames(stairway_geno), , drop = FALSE])
  )

  keep_stairway_individuals <- !is.na(stairway_population)
  stairway_geno <- stairway_geno[keep_stairway_individuals, , drop = FALSE]
  stairway_population <- factor(stairway_population[keep_stairway_individuals])

  # Estimate callable sequence length as the number of RAD tags times 45 bp.
  stairway_RAD_locus <- sub(":.*", "", colnames(stairway_geno))
  stairway_L <- length(unique(stairway_RAD_locus)) * STAIRWAY_RAD_TAG_LENGTH

  stairway_summary <- data.frame()

  for (stairway_pop in levels(stairway_population)) {

    pop_geno <- stairway_geno[
      stairway_population == stairway_pop, , drop = FALSE
    ]
    observed_chromosomes <- 2 * colSums(!is.na(pop_geno))

    # Project onto the even chromosome count covered by approximately 95% of SNPs.
    # Hypergeometric probabilities account for variable missingness across sites.
    nseq <- 2 * floor(as.numeric(quantile(
      observed_chromosomes, .05, na.rm = TRUE
    )) / 2)
    stopifnot(nseq >= 4)

    folded_SFS <- numeric(nseq / 2)
    usable_sites <- 0

    for (locus in seq_len(ncol(pop_geno))) {
      locus_genotypes <- pop_geno[, locus]
      locus_genotypes <- locus_genotypes[!is.na(locus_genotypes)]
      AN <- 2 * length(locus_genotypes)
      AC <- sum(locus_genotypes)

      if (AN < nseq || AC == 0 || AC == AN) next
      usable_sites <- usable_sites + 1

      projected_AC <- 0:nseq
      projected_probability <- dhyper(
        projected_AC,
        AC,
        AN - AC,
        nseq
      )
      folded_class <- pmin(projected_AC, nseq - projected_AC)

      for (frequency_class in seq_len(nseq / 2)) {
        folded_SFS[frequency_class] <-
          folded_SFS[frequency_class] +
          sum(projected_probability[folded_class == frequency_class])
      }
    }

    nrand <- unique(pmax(1, round(
      c(.25, .50, .75, 1) * (nseq / 2 - 1)
    )))
    population_id <- gsub("[^A-Za-z0-9_-]", "_", stairway_pop)
    population_dir <- file.path(stairway_out, population_id)
    dir.create(population_dir, recursive = TRUE, showWarnings = FALSE)
    blueprint_file <- file.path(
      population_dir, paste0(population_id, "_fold.blueprint")
    )

    writeLines(
      c(
        paste0("popid: ", population_id),
        paste0("nseq: ", nseq),
        paste0("L: ", stairway_L),
        "whether_folded: true",
        paste0("SFS: ", paste(format(folded_SFS, scientific = FALSE), collapse = "\t")),
        paste0("pct_training: ", STAIRWAY_TRAINING_PROPORTION),
        paste0("nrand: ", paste(nrand, collapse = "\t")),
        paste0("project_dir: ", population_id),
        paste0(
          "stairway_plot_dir: ",
          file.path(ameiva, "Ameiva_genome", "stairway_plot_v2.1.2", "stairway_plot_es")
        ),
        paste0("ninput: ", STAIRWAY_BOOTSTRAPS),
        paste0("random_seed: ", 2027),
        paste0("mu: ", format(STAIRWAY_MUTATION_RATE, scientific = TRUE)),
        paste0("year_per_generation: ", STAIRWAY_GENERATION_YEARS),
        paste0("plot_title: Ameiva_ameiva_", population_id),
        "xrange: 0,0",
        "yrange: 0,0",
        "xspacing: 2",
        "yspacing: 2",
        "fontsize: 12"
      ),
      blueprint_file
    )

    write.csv(
      data.frame(
        folded_frequency = seq_along(folded_SFS),
        expected_SNP_count = folded_SFS
      ),
      file.path(population_dir, paste0(population_id, "_folded_SFS.csv")),
      row.names = FALSE
    )

    stairway_summary <- rbind(
      stairway_summary,
      data.frame(
        population = stairway_pop,
        individuals = nrow(pop_geno),
        projected_chromosomes = nseq,
        usable_variable_sites = usable_sites,
        RAD_loci = length(unique(stairway_RAD_locus)),
        assumed_callable_bases = stairway_L,
        mutation_rate = STAIRWAY_MUTATION_RATE,
        generation_years = STAIRWAY_GENERATION_YEARS,
        blueprint = blueprint_file
      )
    )
  }

  write.csv(
    stairway_summary,
    file.path(stairway_out, "05_StairwayPlot2_input_audit.csv"),
    row.names = FALSE
  )

  writeLines(
    c(
      "Stairway Plot 2 foi preparado com SFS dobrada para K1-K4.",
      "Os grupos sao o componente sNMF dominante; admistura continua deve ser lembrada.",
      "Nao foi usado o painel MAC>=3 de 14.134 SNPs.",
      paste("Taxa de mutacao assumida:", STAIRWAY_MUTATION_RATE),
      paste("Tempo de geracao assumido:", STAIRWAY_GENERATION_YEARS, "anos"),
      paste("Comprimento de tag assumido:", STAIRWAY_RAD_TAG_LENGTH, "bp"),
      "Repetir com valores alternativos de mu e geracao antes da interpretacao absoluta.",
      "Executar cada blueprint com Stairbuilder e depois o script gerado."
    ),
    file.path(stairway_out, "05_StairwayPlot2_README.txt")
  )
}


# 10.1 Screen environmental predictors for redundancy ----

variables_to_remove <- c(
  # Quarterly variables can change seasonal identity between neighbouring cells
  # (Escobar et al. 2014; Booth 2022). Exclude them before predictor screening.
  "X.CHELSA_bio8_1981.2010_SA",
  "X.CHELSA_bio9_1981.2010_SA",
  "X.CHELSA_bio10_1981.2010_SA",
  "X.CHELSA_bio11_1981.2010_SA",
  "X.CHELSA_bio16_1981.2010_SA",
  "X.CHELSA_bio17_1981.2010_SA",
  "X.CHELSA_bio18_1981.2010_SA",
  "X.CHELSA_bio19_1981.2010_SA",
  "X.CHELSA_npp_1981.2010_SA"
)

stopifnot(all(variables_to_remove %in% names(Variables)))

environment_columns <- grep(
  "^X\\.CHELSA_bio[0-9]+_1981\\.2010_SA$",
  names(Variables),
  value = TRUE
)

environment_columns <- setdiff(environment_columns, variables_to_remove)
Environment_all <- Variables[, environment_columns, drop = FALSE]

names(Environment_all) <- sub(
  "^X\\.CHELSA_(bio[0-9]+)_1981\\.2010_SA$",
  "\\1",
  names(Environment_all)
)

quarter_variables <- paste0("bio", c(8:11, 16:19))
stopifnot(identical(sort(names(Environment_all)), sort(c(paste0("bio", 1:7), paste0("bio", 12:15)))))
stopifnot(!any(quarter_variables %in% names(Environment_all)))
stopifnot(!anyNA(Environment_all))

# Boruta screens climatic predictors of genome-wide PC1.
boruta_data <- data.frame(
  genomic_PC1 = PCA_scores$PC1,
  Environment_all
)

set.seed(2027)
boruta_result <- Boruta(
  genomic_PC1 ~ .,
  data = boruta_data,
  doTrace = 2,
  pValue = .05,
  maxRuns = 1000
)

boruta_final <- TentativeRoughFix(boruta_result)
variables_boruta <- getSelectedAttributes(
  boruta_final,
  withTentative = FALSE
)

variables_boruta
saveRDS(boruta_result, file.path(out, "06_Boruta_environment.rds"))

# Reduce collinearity among predictors retained by Boruta (VIF threshold = 5).
vif_second <- vifstep(Environment_all[, variables_boruta, drop = FALSE], th = 5)
variables_after_vif <- vif_second@results$Variables
variables_after_vif

environment_scaled <- scale(Environment_all[, variables_after_vif, drop = FALSE])
environment_center <- attr(environment_scaled, "scaled:center")
environment_scale <- attr(environment_scaled, "scaled:scale")
Environment_selection <- as.data.frame(environment_scaled)

geography <- data.frame(
  lon = as.numeric(scale(Variables$long)),
  lat = as.numeric(scale(Variables$lat)),
  elevation = as.numeric(scale(Variables$elev.elevation))
)

PCA_condition <- PCA_scores[, c("PC1", "PC2", "PC3")]

# Include K - 1 ancestry coefficients because the K proportions sum to one.
Q_condition <- as.data.frame(Q[, seq_len(bestK - 1), drop = FALSE])
names(Q_condition) <- paste0("Q", seq_len(ncol(Q_condition)))


# 10.2 Select climatic predictors by forward RDA ----

# Follow Josué's workflow: select predictors in an unconditioned RDA,
# then include genomic covariates in the association models.
RDA0 <- rda(geno_complete ~ 1, data = Environment_selection)
RDAfull <- rda(geno_complete ~ ., data = Environment_selection)

forward_josue_file <- file.path(out, paste0(
  "06_forward_RDA_Josue_K", bestK, "_", ENVIRONMENT_VERSION, ".rds"
))
if (USE_SAVED_RDS && file.exists(forward_josue_file)) {
  RDA_forward_josue <- readRDS(forward_josue_file)
} else {
  set.seed(2027)
  RDA_forward_josue <- ordiR2step(
    RDA0,
    scope = formula(RDAfull),
    direction = "forward",
    Pin = .01,
    R2permutations = 1000,
    R2scope = TRUE,
    permutations = 999,
    trace = TRUE
  )
  saveRDS(RDA_forward_josue, forward_josue_file)
}

variables_final <- attr(terms(RDA_forward_josue), "term.labels")
stopifnot(length(variables_final) > 0)
Environment <- Environment_selection[, variables_final, drop = FALSE]

write.csv(
  data.frame(
    method = "Boruta -> VIF -> ordiR2step unconditioned (Josue)",
    selected = paste(variables_final, collapse = ";")
  ),
  file.path(out, "06_environment_forward_selection.csv"),
  row.names = FALSE
)


# 11. Fit partial RDA models and compare conditioning schemes ----

RDA_data <- cbind(Environment, geography, PCA_condition, Q_condition)
snmf_terms <- paste(names(Q_condition), collapse = " + ")

# The response contains individuals by SNPs; climate is the explanatory matrix.
# Begin with the model conditioned on genomic PC1–PC3.
formula_josue <- as.formula(paste(
  "geno_complete ~", paste(names(Environment), collapse = " + "),
  "+ Condition(PC1 + PC2 + PC3)"
))
RDA_josue <- rda(formula_josue, data = RDA_data)

# Compare unconditioned RDA with two partial RDA models that account for
# geography and either genome-wide PCs or sNMF ancestry coefficients.
RDA_simple <- rda(geno_complete ~ ., data = Environment)
formula_geo <- as.formula(paste(
  "geno_complete ~", paste(names(Environment), collapse = " + "),
  "+ Condition(PC1 + PC2 + PC3 + lon + lat)"
))
RDA_geo <- rda(formula_geo, data = RDA_data)
formula_snmf_geo <- as.formula(paste(
  "geno_complete ~", paste(names(Environment), collapse = " + "),
  "+ Condition(", snmf_terms, "+ lon + lat)"
))
RDA_snmf_geo <- rda(formula_snmf_geo, data = RDA_data)

# Elevation may capture isolation and climatic variation; include it in sensitivity models.
formula_geo_elevation <- as.formula(paste(
  "geno_complete ~", paste(names(Environment), collapse = " + "),
  "+ Condition(PC1 + PC2 + PC3 + lon + lat + elevation)"
))
RDA_geo_elevation <- rda(formula_geo_elevation, data = RDA_data)
formula_snmf_geo_elevation <- as.formula(paste(
  "geno_complete ~", paste(names(Environment), collapse = " + "),
  "+ Condition(", snmf_terms, "+ lon + lat + elevation)"
))
RDA_snmf_geo_elevation <- rda(formula_snmf_geo_elevation, data = RDA_data)

RDA_models <- list(
  Josue_PC1_PC3 = RDA_josue,
  PCA_geo = RDA_geo,
  sNMF_geo = RDA_snmf_geo,
  PCA_geo_elevation_sensitivity = RDA_geo_elevation,
  sNMF_geo_elevation_sensitivity = RDA_snmf_geo_elevation,
  simple_sensitivity = RDA_simple
)

# Confirm the response retains the 245 individuals and 14,134 SNPs.
stopifnot(all(vapply(
  RDA_models,
  function(model) identical(dim(model$Ybar), dim(geno_complete)),
  logical(1)
)))

RDA_design_audit <- data.frame(
  model = names(RDA_models),
  response = "allelic dosage matrix",
  response_rows_individuals = nrow(geno_complete),
  response_columns_SNPs = ncol(geno_complete),
  explanatory_matrix = paste(length(variables_final), "standardized climate predictors"),
  conditioned_covariates = c(
    "genomic PC1-PC3",
    "genomic PC1-PC3 + longitude + latitude",
    paste0("sNMF K=", bestK, " Q1-Q3 + longitude + latitude"),
    "genomic PC1-PC3 + longitude + latitude + elevation; sensitivity",
    paste0("sNMF K=", bestK, " Q1-Q3 + longitude + latitude + elevation; sensitivity"),
    "none; confounding sensitivity only"
  ),
  formula = vapply(
    RDA_models,
    function(model) paste(deparse(formula(model)), collapse = " "),
    character(1)
  )
)
write.csv(
  RDA_design_audit,
  file.path(out, "07_RDA_response_predictor_audit.csv"),
  row.names = FALSE
)

for (model_name in names(RDA_models)) {
  model <- RDA_models[[model_name]]
  model_file <- file.path(out, paste0("07_RDA_", model_name, ".rds"))
  saveRDS(model, model_file)

  print(model_name)
  print(RsquareAdj(model))
  print(vif.cca(model))

  # Use 999 permutations for the global test and 199 for optional axis and margin diagnostics.
  global_test_file <- file.path(
    out,
    paste0("07_RDA_", model_name, "_K", bestK, "_", ENVIRONMENT_VERSION, "_global_test.rds")
  )
  test_global <- if (USE_SAVED_RDS && file.exists(global_test_file)) {
    readRDS(global_test_file)
  } else {
    anova(model, permutations = 999)
  }
  print(test_global)
  saveRDS(test_global, global_test_file)
  rm(test_global)
  gc()

  if (RUN_PARTITIONED_RDA_TESTS) {
    # Axis and marginal tests have substantial memory requirements for this SNP panel.
    test_axis <- anova(model, by = "axis", permutations = 199)
    print(test_axis)
    saveRDS(test_axis, file.path(out, paste0("07_RDA_", model_name, "_axis_test.rds")))
    rm(test_axis)
    gc()

    test_margin <- anova(model, by = "margin", permutations = 199)
    print(test_margin)
    saveRDS(test_margin, file.path(out, paste0("07_RDA_", model_name, "_margin_test.rds")))
    rm(test_margin)
    gc()
  }
}

# Plot individual scores and climatic predictor vectors in partial RDA space.
RDA_geo_sites <- as.data.frame(scores(
  RDA_geo, display = "sites", choices = 1:2, scaling = 2
))
RDA_geo_arrows <- as.data.frame(scores(
  RDA_geo, display = "bp", choices = 1:2, scaling = 2
))
names(RDA_geo_sites)[1:2] <- c("RDA1", "RDA2")
names(RDA_geo_arrows)[1:2] <- c("RDA1", "RDA2")
RDA_geo_sites$ancestry_cluster <- factor(max.col(Q))
RDA_geo_arrows$environment <- rownames(RDA_geo_arrows)

RDA_arrow_multiplier <- .70 * min(
  max(abs(RDA_geo_sites$RDA1), na.rm = TRUE) /
    max(abs(RDA_geo_arrows$RDA1), na.rm = TRUE),
  max(abs(RDA_geo_sites$RDA2), na.rm = TRUE) /
    max(abs(RDA_geo_arrows$RDA2), na.rm = TRUE)
)
RDA_geo_arrows$RDA1 <- RDA_geo_arrows$RDA1 * RDA_arrow_multiplier
RDA_geo_arrows$RDA2 <- RDA_geo_arrows$RDA2 * RDA_arrow_multiplier
RDA_axis_fraction <- RDA_geo$CCA$eig / sum(RDA_geo$CCA$eig)

plot_partial_RDA <- ggplot(
  RDA_geo_sites, aes(RDA1, RDA2, colour = ancestry_cluster)
) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = .35) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = .35) +
  geom_point(size = 2.2, alpha = .82) +
  geom_segment(
    data = RDA_geo_arrows,
    aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(2.2, "mm")),
    colour = idv_orange_dark, linewidth = .75
  ) +
  geom_text(
    data = RDA_geo_arrows,
    aes(RDA1, RDA2, label = environment),
    inherit.aes = FALSE,
    colour = idv_orange_dark, fontface = "bold", size = 3,
    nudge_y = .025 * diff(range(RDA_geo_sites$RDA2))
  ) +
  scale_colour_manual(values = idv_discrete(bestK)) +
  coord_equal() +
  theme_ameiva() +
  labs(
    title = "Structure-controlled genotype-environment association",
    subtitle = "Response: 245 x 14,134 allelic matrix; arrows: climatic predictors",
    x = paste0("RDA1 (", round(100 * RDA_axis_fraction[1], 1), "% constrained)"),
    y = paste0("RDA2 (", round(100 * RDA_axis_fraction[2], 1), "% constrained)"),
    colour = "sNMF cluster"
  )

ggsave(
  file.path(figures, "31_partial_RDA_genotype_response.png"),
  plot_partial_RDA, width = 8, height = 6,
  dpi = 400, bg = "white"
)


# 12. Detect partial RDA candidates using robust locus-loading distances ----

rdadapt <- function(rda_object, K) {
  loadings <- rda_object$CCA$v[, seq_len(K), drop = FALSE]
  locus_ids <- rownames(loadings)
  if (is.null(locus_ids)) locus_ids <- colnames(geno_complete)
  stopifnot(length(locus_ids) == nrow(loadings))

  # Store SNP identifiers before scaling loadings to retain the link to the VCF.
  loadings <- scale(loadings)
  max_abs_loading_z <- apply(abs(loadings), 1, max)
  robust_distance <- covRob(
    loadings,
    distance = TRUE,
    na.action = na.omit,
    estim = "pairwiseGK"
  )$dist
  lambda <- median(robust_distance) / qchisq(.5, df = K)
  p_value <- pchisq(robust_distance / lambda, K, lower.tail = FALSE)
  data.frame(
    stacks_id = locus_ids,
    p_value = p_value,
    q_value = qvalue::qvalue(p_value)$qvalues,
    robust_distance = robust_distance,
    max_abs_loading_z = max_abs_loading_z
  )
}

K_axes <- min(3, RDA_josue$CCA$rank)
rda_scans <- lapply(RDA_models, function(x) rdadapt(x, min(K_axes, x$CCA$rank)))

for (model_name in names(rda_scans)) {
  scan <- rda_scans[[model_name]]
  # Use Bonferroni-adjusted significance for candidate selection; retain q-values for comparison.
  bonferroni_threshold <- .01 / nrow(scan)
  scan$bonferroni_threshold <- bonferroni_threshold
  scan$candidate_Bonferroni01 <- scan$p_value < bonferroni_threshold
  scan$candidate_3SD <- scan$max_abs_loading_z >= 3
  rda_scans[[model_name]] <- scan
  write.csv(scan, file.path(out, paste0("08_RDA_locus_scan_", model_name, ".csv")), row.names = FALSE)
}

candidates_RDA <- rda_scans$Josue_PC1_PC3$stacks_id[
  rda_scans$Josue_PC1_PC3$candidate_Bonferroni01
]
candidates_RDA_pca_geo <- rda_scans$PCA_geo$stacks_id[
  rda_scans$PCA_geo$candidate_Bonferroni01
]
candidates_RDA_snmf_geo <- rda_scans$sNMF_geo$stacks_id[
  rda_scans$sNMF_geo$candidate_Bonferroni01
]
candidates_RDA_material <- intersect(
  candidates_RDA_pca_geo,
  candidates_RDA_snmf_geo
)
candidates_RDA_material_strict <- intersect(
  rda_scans$PCA_geo$stacks_id[
    rda_scans$PCA_geo$candidate_Bonferroni01 &
      rda_scans$PCA_geo$candidate_3SD
  ],
  rda_scans$sNMF_geo$stacks_id[
    rda_scans$sNMF_geo$candidate_Bonferroni01 &
      rda_scans$sNMF_geo$candidate_3SD
  ]
)
write.csv(data.frame(stacks_id = candidates_RDA), file.path(out, "08_RDA_candidates_Josue.csv"), row.names = FALSE)
write.csv(data.frame(stacks_id = candidates_RDA_material), file.path(out, "08_RDA_candidates_landscape_consensus.csv"), row.names = FALSE)
write.csv(data.frame(stacks_id = candidates_RDA_material_strict), file.path(out, "08_RDA_candidates_landscape_strict_3SD.csv"), row.names = FALSE)

# Assess linkage sensitivity by retaining the lowest-P SNP per RAD locus.
rad_id_from_stacks <- function(x) sub(":.*$", "", x)
best_per_rad <- function(scan, candidate_column = "candidate_Bonferroni01") {
  selected <- scan[scan[[candidate_column]], , drop = FALSE]
  if (!nrow(selected)) return(selected)
  selected$rad_locus <- rad_id_from_stacks(selected$stacks_id)
  selected <- selected[order(selected$rad_locus, selected$p_value), , drop = FALSE]
  selected[!duplicated(selected$rad_locus), , drop = FALSE]
}
write.csv(
  best_per_rad(rda_scans$PCA_geo),
  file.path(out, "08_RDA_candidates_PCA_geo_Bonferroni_best_per_RAD.csv"),
  row.names = FALSE
)
write.csv(
  best_per_rad(rda_scans$sNMF_geo),
  file.path(out, "08_RDA_candidates_sNMF_geo_Bonferroni_best_per_RAD.csv"),
  row.names = FALSE
)


# 13. Test each climatic predictor with LFMM2 ----

if (RUN_LFMM) {

  # Test each selected climatic predictor against the same sNMF-imputed allele matrix.
  stopifnot(!anyNA(geno_complete))
  stopifnot(identical(rownames(geno_complete), rownames(Environment)))
  stopifnot(identical(colnames(geno_complete), colnames(geno)))

  # Set the latent-factor count to the selected sNMF K. Optional fits at K = 2–8
  # compare genomic inflation factors without selecting K by association P-values.
  K_LFMM2 <- bestK
  K_sensitivity <- 2:8
  K_sensitivity <- K_sensitivity[K_sensitivity < nrow(geno_complete)]
  RUN_LFMM2_K_SENSITIVITY <- FALSE
  stopifnot(K_LFMM2 %in% K_sensitivity)

  lfmm2_results <- vector("list", ncol(Environment))
  lfmm2_gif_sensitivity <- vector("list", ncol(Environment))

  for (j in seq_len(ncol(Environment))) {
    predictor_name <- colnames(Environment)[j]
    X_predictor <- as.matrix(Environment[, j, drop = FALSE])

    # Fit the prespecified latent-factor model and calibrate the tests.
    set.seed(2027 + j)
    fit_LFMM2 <- LEA::lfmm2(
      input = geno_complete,
      env = X_predictor,
      K = K_LFMM2
    )
    test_LFMM2 <- LEA::lfmm2.test(
      object = fit_LFMM2,
      input = geno_complete,
      env = X_predictor,
      full = FALSE,
      genomic.control = TRUE
    )

    p_LFMM2 <- as.numeric(test_LFMM2$pvalues)
    lfmm2_results[[j]] <- data.frame(
      stacks_id = colnames(geno_complete),
      environmental_predictor = predictor_name,
      K_sNMF = K_LFMM2,
      GIF = as.numeric(test_LFMM2$gif[1]),
      z_score = as.numeric(test_LFMM2$zscores),
      p_value = p_LFMM2,
      q_value = p.adjust(p_LFMM2, method = "BH")
    )

    # Compare genomic inflation factors across the requested K values.
    K_to_report <- if (RUN_LFMM2_K_SENSITIVITY) K_sensitivity else K_LFMM2
    gif_values <- numeric(length(K_to_report))
    for (k_index in seq_along(K_to_report)) {
      set.seed(202700 + 100 * j + K_to_report[k_index])
      fit_sensitivity <- LEA::lfmm2(
        input = geno_complete,
        env = X_predictor,
        K = K_to_report[k_index]
      )
      test_sensitivity <- LEA::lfmm2.test(
        object = fit_sensitivity,
        input = geno_complete,
        env = X_predictor,
        full = FALSE,
        genomic.control = TRUE
      )
      gif_values[k_index] <- as.numeric(test_sensitivity$gif[1])
      rm(fit_sensitivity, test_sensitivity)
      gc(verbose = FALSE)
    }
    lfmm2_gif_sensitivity[[j]] <- data.frame(
      environmental_predictor = predictor_name,
      K_sNMF = K_to_report,
      GIF = gif_values,
      analysis_K = K_to_report == K_LFMM2
    )

    rm(fit_LFMM2, test_LFMM2)
    gc(verbose = FALSE)
  }

  lfmm2_results <- do.call(rbind, lfmm2_results)
  lfmm2_gif_sensitivity <- do.call(rbind, lfmm2_gif_sensitivity)
  lfmm2_candidates <- lfmm2_results[lfmm2_results$q_value <= .05, , drop = FALSE]

  write.csv(lfmm2_results, file.path(out, "09_LFMM2_all_14134.csv"), row.names = FALSE)
  write.csv(lfmm2_candidates, file.path(out, "09_LFMM2_candidates_FDR05.csv"), row.names = FALSE)
  write.csv(lfmm2_gif_sensitivity, file.path(out, "09_LFMM2_GIF_K_sensitivity.csv"), row.names = FALSE)

  candidates_LFMM <- unique(lfmm2_candidates$stacks_id)
  # The shared linear candidate set is the intersection of the landscape
  # partial RDA consensus and LFMM2 candidates.
  linear_adaptive_loci <- intersect(
    candidates_RDA_material,
    candidates_LFMM
  )
  candidates_consensus <- linear_adaptive_loci
  candidates_evidence_union <- union(candidates_RDA_material, candidates_LFMM)

  evidence <- data.frame(stacks_id = colnames(geno_complete))
  evidence$RDA_Josue <- evidence$stacks_id %in% candidates_RDA
  evidence$RDA_landscape_consensus <- evidence$stacks_id %in% candidates_RDA_material
  for (j in seq_len(ncol(Environment))) {
    predictor_candidates <- lfmm2_candidates$stacks_id[
      lfmm2_candidates$environmental_predictor == colnames(Environment)[j]
    ]
    evidence[[paste0("LFMM2_", colnames(Environment)[j])]] <-
      evidence$stacks_id %in% predictor_candidates
  }
  evidence$evidence_count <- rowSums(evidence[, -1])
  evidence$tier <- cut(
    evidence$evidence_count,
    breaks = c(-1, 0, 1, 2, Inf),
    labels = c("background", "single_method", "two_methods", "three_or_more_methods")
  )

  write.csv(evidence, file.path(out, "09_candidate_evidence_tiers.csv"), row.names = FALSE)
  write.csv(
    evidence[evidence$stacks_id %in% linear_adaptive_loci, , drop = FALSE],
    file.path(out, "09_linear_adaptive_loci_RDA_and_LFMM.csv"),
    row.names = FALSE
  )
}


# 14. Annotate SNPs and compare functional classes ----

# Classify SNPs across the complete 14,134-marker panel to define the background
# for functional enrichment. Use annotations transferred to the Ameiva assembly.
annotation_14134_file <- file.path(out, "10_annotation_14134.csv")

if (!USE_SAVED_RDS || !file.exists(annotation_14134_file)) {

  transferred_gff <- file.path(
    hd,
    "work",
    "Aspidoscelis_on_Ameiva_GCA_046270845.1.gff3"
  )
  stopifnot(file.exists(transferred_gff))

  gff <- read.delim(
    transferred_gff,
    comment.char = "#",
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    col.names = c(
      "seqid", "source", "type", "start", "end",
      "score", "strand", "phase", "attributes"
    )
  )

  gff_gr <- GRanges(
    gff$seqid,
    IRanges(gff$start, gff$end),
    strand = gff$strand,
    type = gff$type,
    attributes = gff$attributes
  )

  loci_gr <- GRanges(
    qc_loci$scaffold,
    IRanges(qc_loci$position, width = 1),
    strand = "*",
    stacks_id = qc_loci$stacks_id
  )

  genes <- gff_gr[tolower(mcols(gff_gr)$type) == "gene"]
  exons <- gff_gr[tolower(mcols(gff_gr)$type) == "exon"]
  cds <- gff_gr[tolower(mcols(gff_gr)$type) == "cds"]

  in_gene <- unique(queryHits(findOverlaps(loci_gr, genes, ignore.strand = TRUE)))
  in_exon <- unique(queryHits(findOverlaps(loci_gr, exons, ignore.strand = TRUE)))
  in_cds <- unique(queryHits(findOverlaps(loci_gr, cds, ignore.strand = TRUE)))

  annotation_14134 <- qc_loci
  annotation_14134$functional_class <- "intergenic"
  annotation_14134$functional_class[in_gene] <- "intronic"
  annotation_14134$functional_class[in_exon] <- "exonic_non_CDS"
  annotation_14134$functional_class[in_cds] <- "CDS"

  gene_hits <- findOverlaps(loci_gr, genes, ignore.strand = TRUE)
  gene_for_locus <- rep(NA_integer_, length(loci_gr))
  gene_for_locus[queryHits(gene_hits)] <- subjectHits(gene_hits)

  without_gene <- which(is.na(gene_for_locus))
  nearest <- distanceToNearest(loci_gr[without_gene], genes, ignore.strand = TRUE)
  gene_for_locus[without_gene[queryHits(nearest)]] <- subjectHits(nearest)

  gene_attributes <- as.character(mcols(genes)$attributes[gene_for_locus])
  gene_id <- sub(".*(?:^|;)ID=([^;]+).*", "\\1", gene_attributes, perl = TRUE)
  gene_name <- sub(
    ".*(?:^|;)(?:Name|gene_name)=([^;]+).*",
    "\\1",
    gene_attributes,
    perl = TRUE
  )

  annotation_14134$gene_id <- ifelse(
    gene_id != gene_attributes,
    gene_id,
    NA_character_
  )
  annotation_14134$gene_name <- ifelse(
    gene_name != gene_attributes,
    gene_name,
    NA_character_
  )

  annotation_14134$distance_to_gene <- NA_integer_
  annotation_14134$distance_to_gene[in_gene] <- 0L
  annotation_14134$distance_to_gene[
    without_gene[queryHits(nearest)]
  ] <- mcols(nearest)$distance

  write.csv(annotation_14134, annotation_14134_file, row.names = FALSE)
}

annotation_14134 <- read.csv(annotation_14134_file, check.names = FALSE)
candidate_annotation <- annotation_14134[
  annotation_14134$stacks_id %in% candidates_evidence_union,
]
linear_adaptive_annotation <- annotation_14134[
  annotation_14134$stacks_id %in% linear_adaptive_loci,
]

write.csv(
  candidate_annotation,
  file.path(out, "10_candidate_annotation.csv"),
  row.names = FALSE
)
write.csv(
  candidate_annotation,
  file.path(out, "10_candidate_annotation_RDA_or_LFMM.csv"),
  row.names = FALSE
)
write.csv(
  linear_adaptive_annotation,
  file.path(out, "10_linear_adaptive_annotation_RDA_and_LFMM.csv"),
  row.names = FALSE
)

annotation_summary <- as.data.frame(table(
  universe = annotation_14134$functional_class
))
candidate_summary <- as.data.frame(table(
  candidates = candidate_annotation$functional_class
))
linear_adaptive_summary <- as.data.frame(table(
  functional_class = linear_adaptive_annotation$functional_class
))

linear_adaptive_enrichment <- do.call(rbind, lapply(
  unique(annotation_14134$functional_class),
  function(functional_class) {
    in_class <- annotation_14134$functional_class == functional_class
    is_adaptive <- annotation_14134$stacks_id %in% linear_adaptive_loci
    contingency <- matrix(c(
      sum(in_class & is_adaptive),
      sum(in_class & !is_adaptive),
      sum(!in_class & is_adaptive),
      sum(!in_class & !is_adaptive)
    ), nrow = 2, byrow = TRUE)
    test <- fisher.test(contingency)
    data.frame(
      functional_class = functional_class,
      adaptive_count = sum(in_class & is_adaptive),
      background_count = sum(in_class & !is_adaptive),
      odds_ratio = unname(test$estimate),
      fisher_p = test$p.value
    )
  }
))
linear_adaptive_enrichment$fdr <- p.adjust(
  linear_adaptive_enrichment$fisher_p,
  method = "BH"
)

write.csv(annotation_summary, file.path(out, "10_annotation_universe_summary.csv"), row.names = FALSE)
write.csv(candidate_summary, file.path(out, "10_annotation_candidate_summary.csv"), row.names = FALSE)
write.csv(
  linear_adaptive_summary,
  file.path(out, "10_linear_adaptive_functional_summary.csv"),
  row.names = FALSE
)
write.csv(
  linear_adaptive_enrichment,
  file.path(out, "10_linear_adaptive_functional_enrichment.csv"),
  row.names = FALSE
)

functional_plot_data <- linear_adaptive_summary
names(functional_plot_data) <- c("functional_class", "candidate_SNPs")
functional_plot_data$functional_class <- factor(
  functional_plot_data$functional_class,
  levels = c("intergenic", "intronic", "exonic_non_CDS", "CDS")
)
plot_candidate_function <- ggplot(
  functional_plot_data,
  aes(functional_class, candidate_SNPs, fill = functional_class)
) +
  geom_col(width = .72) +
  geom_text(
    aes(label = candidate_SNPs),
    vjust = -.35, colour = idv_purple_dark, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    intergenic = idv_purple_dark,
    intronic = idv_purple_light,
    exonic_non_CDS = idv_orange_light,
    CDS = idv_orange_dark
  ), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  theme_ameiva() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(
    title = "Functional context of linear candidates",
    subtitle = paste(length(candidates_consensus), "SNPs supported jointly by landscape RDA and LFMM"),
    x = NULL, y = "Candidate SNPs"
  )

ggsave(
  file.path(figures, "34_candidate_functional_context.png"),
  plot_candidate_function, width = 7, height = 5.5,
  dpi = 400, bg = "white"
)


# 15. Explore Gradient Forest turnover for shared linear candidates ----

# This candidate-panel GF uses the four-variable linear association workflow.
# The separate full-panel GF script models all 14,134 SNPs.

GF_loci <- candidates_consensus
GF_loci <- intersect(GF_loci, colnames(geno_complete))

if (RUN_GRADIENT_FOREST) {

  stopifnot(length(GF_loci) >= 2)

  GF_response <- geno_complete[, GF_loci, drop = FALSE]
  GF_response_names <- paste0("SNP", seq_len(ncol(GF_response)))
  GF_map <- data.frame(
    response_name = GF_response_names,
    stacks_id = colnames(GF_response)
  )
  colnames(GF_response) <- GF_response_names

  GF_data <- data.frame(Environment, GF_response, check.names = FALSE)
  GF_model_file <- file.path(out, paste0("11_Gradient_Forest_", ENVIRONMENT_VERSION, ".rds"))
  GF_map_file <- file.path(out, "11_Gradient_Forest_SNP_map.csv")
  reuse_GF <- USE_SAVED_RDS && file.exists(GF_model_file) && file.exists(GF_map_file) &&
    identical(read.csv(GF_map_file)$stacks_id, GF_map$stacks_id)

  if (reuse_GF) {
    gradient_forest <- readRDS(GF_model_file)
  } else {
    set.seed(2027)
    gradient_forest <- gradientForest(
      data = GF_data,
      predictor.vars = names(Environment),
      response.vars = GF_response_names,
      ntree = 1000,
      corr.threshold = .50,
      compact = TRUE,
      trace = TRUE,
      check.names = FALSE
    )
    saveRDS(gradient_forest, GF_model_file)
  }

  write.csv(GF_map, file.path(out, "11_Gradient_Forest_SNP_map.csv"), row.names = FALSE)

  GF_importance <- data.frame(
    environment = names(gradient_forest$overall.imp),
    importance = as.numeric(gradient_forest$overall.imp)
  )
  GF_importance <- GF_importance[order(GF_importance$importance, decreasing = TRUE),]
  write.csv(GF_importance, file.path(out, "11_Gradient_Forest_importance.csv"), row.names = FALSE)

  plot_GF_importance <- ggplot(
    GF_importance,
    aes(reorder(environment, importance), importance, fill = importance)
  ) +
    geom_col(width = .72) +
    scale_fill_gradientn(colours = idv_sequential, guide = "none") +
    coord_flip() +
    theme_ameiva() +
    labs(
      title = "Environmental importance in the adaptive mosaic",
      subtitle = paste(length(GF_loci), "loci in landscape-RDA x LFMM consensus"),
      x = NULL,
      y = "Gradient Forest importance"
    )

  ggsave(file.path(figures, "03_Gradient_Forest_importance.png"), plot_GF_importance,
         width = 8, height = 6, dpi = 400, bg = "white")
}


# 16. Model genome-wide nonlinear associations with RDAforest ----

# Follow the wolf vignette in RDAforest 2.6.10 (https://github.com/z0on/RDA-forest).
# Use 1 minus genotype correlation as genomic dissimilarity, with geographic
# axes and genomic PC1–PC3 as covariates.

genomic_distance <- 1 - cor(t(geno_complete), use = "pairwise.complete.obs")
genomic_distance <- (genomic_distance + t(genomic_distance)) / 2
diag(genomic_distance) <- 0

gcd_result <- RDAforest::gcd.dist(Variables[, c("long", "lat")])
gcd <- gcd_result$gcd.coords
geographic_distance <- as.matrix(gcd_result$gcd.distances)

set.seed(2027)
RDAforest_IBD_test <- protest(
  capscale(geographic_distance ~ 1),
  capscale(genomic_distance ~ 1),
  permutations = 999
)

write.csv(
  data.frame(
    protest_correlation = RDAforest_IBD_test$t0,
    permutation_P = RDAforest_IBD_test$signif,
    permutations = 999
  ),
  file.path(out, "12_RDAforest_IBD_test.csv"),
  row.names = FALSE
)

RDAforest_covariates <- data.frame(
  geographic_axis1 = gcd[, 1],
  geographic_axis2 = gcd[, 2],
  PC1 = PCA_condition$PC1,
  PC2 = PCA_condition$PC2,
  PC3 = PCA_condition$PC3
)

RDAforest_covariates <- as.data.frame(scale(RDAforest_covariates))
RDAforest_environment <- as.data.frame(scale(Environment))

RDAforest_ordination <- capscale(
  genomic_distance ~ 1 + Condition(as.matrix(RDAforest_covariates))
)
RDAforest_max_pcs <- min(40, length(RDAforest_ordination$CA$eig))
RDAforest_residual_scores <- scores(
  RDAforest_ordination,
  display = "sites",
  choices = seq_len(RDAforest_max_pcs),
  scaling = 1
)
RDAforest_axis_variance <- apply(RDAforest_residual_scores, 2, var)
RDAforest_candidate_elbows <- 2:min(25, RDAforest_max_pcs - 1)
RDAforest_piecewise_sse <- vapply(
  RDAforest_candidate_elbows,
  function(k) {
    first <- data.frame(
      axis = seq_len(k),
      variance = RDAforest_axis_variance[seq_len(k)]
    )
    second_index <- k:RDAforest_max_pcs
    second <- data.frame(
      axis = second_index,
      variance = RDAforest_axis_variance[second_index]
    )
    sum(residuals(lm(log10(variance) ~ axis, data = first))^2) +
      sum(residuals(lm(log10(variance) ~ axis, data = second))^2)
  },
  numeric(1)
)
RDAforest_top_pcs <- RDAforest_candidate_elbows[
  which.min(RDAforest_piecewise_sse)
]
RDAforest_elbow_pcs <- RDAforest_top_pcs
stopifnot(RDAforest_elbow_pcs == 8)
if (is.finite(RDAFOREST_PCS_OVERRIDE) && RDAFOREST_PCS_OVERRIDE > 0) {
  RDAforest_top_pcs <- RDAFOREST_PCS_OVERRIDE
}
RDAforest_axis_audit <- data.frame(
  axis = seq_len(RDAforest_max_pcs),
  variance = RDAforest_axis_variance,
  percent_among_first40 = 100 * RDAforest_axis_variance /
    sum(RDAforest_axis_variance),
  cumulative_percent_among_first40 = 100 * cumsum(RDAforest_axis_variance) /
    sum(RDAforest_axis_variance),
  retained_by_piecewise_elbow = seq_len(RDAforest_max_pcs) <=
    RDAforest_top_pcs
)
write.csv(
  RDAforest_axis_audit,
  file.path(out, "12_RDAforest_axis_variance_elbow.csv"),
  row.names = FALSE
)
RDAforest_residual_distance <- as.matrix(dist(RDAforest_residual_scores))
set.seed(2029)
RDAforest_residual_IBD_test <- mantel(
  as.dist(RDAforest_residual_distance),
  as.dist(geographic_distance),
  method = "spearman",
  permutations = 999
)
write.csv(
  data.frame(
    diagnostic = c("raw genomic distance", "conditioned residual ordination"),
    method = c("PROTEST", "Spearman Mantel"),
    statistic = c(
      RDAforest_IBD_test$t0,
      unname(RDAforest_residual_IBD_test$statistic)
    ),
    permutation_P = c(
      RDAforest_IBD_test$signif,
      RDAforest_residual_IBD_test$signif
    ),
    permutations = 999
  ),
  file.path(out, "12_RDAforest_raw_and_residual_IBD_test.csv"),
  row.names = FALSE
)

if (RUN_RDAFOREST) {

  exploratory_file <- file.path(out, paste0(
    "12_RDAforest_exploratory_v2610_", ENVIRONMENT_VERSION, ".rds"
  ))
  mtry_file <- file.path(out, paste0(
    "12_RDAforest_mtry_selection_v2610_", ENVIRONMENT_VERSION,
    "_pcs", RDAforest_top_pcs, ".rds"
  ))
  rdaforest_file <- file.path(out, paste0(
    "12_RDAforest_fit_v2610_", ENVIRONMENT_VERSION,
    "_pcs", RDAforest_top_pcs, ".rds"
  ))

  if (USE_SAVED_RDS && file.exists(exploratory_file)) {
    RDAforest_exploratory <- readRDS(exploratory_file)
  } else {
    set.seed(2027)
    RDAforest_exploratory <- RDAforest::makeGF(
      RDAforest_ordination,
      RDAforest_environment,
      pcs2keep = seq_len(RDAforest_max_pcs),
      ntrees = 1500
    )
    saveRDS(RDAforest_exploratory, exploratory_file)
  }

  RDAforest_PC_prediction <- data.frame(
    PC = names(RDAforest_exploratory$result),
    OOB_R2 = as.numeric(RDAforest_exploratory$result)
  )
  write.csv(
    RDAforest_PC_prediction,
    file.path(out, "12_RDAforest_PC_predictability.csv"),
    row.names = FALSE
  )

  if (USE_SAVED_RDS && file.exists(mtry_file) && file.exists(rdaforest_file)) {
    mtry_selection <- readRDS(mtry_file)
    rdaforest_fit <- readRDS(rdaforest_file)
    selected_RDAforest <- mtry_selection$goodvars
  } else {
    set.seed(2027)
    mtry_selection <- RDAforest::mtrySelJack(
      Y = genomic_distance,
      X = RDAforest_environment,
      covariates = RDAforest_covariates,
      nreps = 31,
      oob = .20,
      prop.positive.cutoff = .50,
      importance.cutoff = .10,
      importance.type = "GF",
      top.pcs = RDAforest_top_pcs,
      ntrees = 1500
    )
    selected_RDAforest <- mtry_selection$goodvars
    stopifnot(length(selected_RDAforest) > 0)
    set.seed(2028)
    rdaforest_fit <- RDAforest::ordinationJackknife(
      Y = genomic_distance,
      X = RDAforest_environment[, selected_RDAforest, drop = FALSE],
      covariates = RDAforest_covariates,
      nreps = 25,
      oob = .20,
      top.pcs = RDAforest_top_pcs,
      ntrees = 1500
    )
    saveRDS(mtry_selection, mtry_file)
    saveRDS(rdaforest_fit, rdaforest_file)
  }

  stopifnot(mtry_selection$failed.reps == 0)
  stopifnot(rdaforest_fit$failed.reps == 0)

  write.csv(
    mtry_selection$delta,
    file.path(out, "12_RDAforest_mtry_delta.csv"),
    row.names = FALSE
  )
  write.csv(
    mtry_selection$prop.positive,
    file.path(out, "12_RDAforest_mtry_proportion_positive.csv"),
    row.names = FALSE
  )
  write.csv(
    rdaforest_fit$all.importances,
    file.path(out, "12_RDAforest_importance.csv"),
    row.names = FALSE
  )

  RDAforest_importance <- aggregate(
    importance ~ variable,
    rdaforest_fit$all.importances,
    median
  )
  RDAforest_importance <- RDAforest_importance[
    order(RDAforest_importance$importance, decreasing = TRUE),
  ]

  plot_RDAforest <- ggplot(
    RDAforest_importance,
    aes(reorder(variable, importance), importance, fill = importance)
  ) +
    geom_col(width = .72) +
    scale_fill_gradientn(colours = idv_sequential, guide = "none") +
    coord_flip() +
    theme_ameiva() +
    labs(
      title = "RDAforest: nonlinear genomic covariance",
      subtitle = "Population structure and geographic position treated as covariates",
      x = NULL,
      y = "Median jackknife importance"
    )

  ggsave(file.path(figures, "04_RDAforest_importance.png"), plot_RDAforest,
         width = 8, height = 6, dpi = 400, bg = "white")

  # Compare geographic and genomic distances before and after conditioning.
  pair_index <- lower.tri(genomic_distance)
  ibd_table <- rbind(
    data.frame(
      geographic_km = geographic_distance[pair_index],
      genomic_distance = genomic_distance[pair_index],
      stage = "Raw genomic distance"
    ),
    data.frame(
      geographic_km = geographic_distance[pair_index],
      genomic_distance = RDAforest_residual_distance[pair_index],
      stage = "After conditioning"
    )
  )
  plot_IBD <- ggplot(ibd_table, aes(geographic_km, genomic_distance)) +
    geom_point(alpha = .055, size = .45, colour = idv_purple) +
    stat_summary_bin(
      fun = median, geom = "line", bins = 40,
      colour = idv_orange, linewidth = 1
    ) +
    facet_wrap(~stage, scales = "free_y") +
    theme_ameiva() +
    labs(
      title = "RDAforest diagnostic: isolation by distance",
      subtitle = paste0(
        "PROTEST r raw = ", round(RDAforest_IBD_test$t0, 3),
        "; residual Mantel rho = ",
        round(unname(RDAforest_residual_IBD_test$statistic), 3)
      ),
      x = "Great-circle distance (km)",
      y = "Genomic dissimilarity (panel-specific units)"
    )
  ggsave(file.path(figures, "14_RDAforest_IBD_diagnostic.png"), plot_IBD,
         width = 8, height = 6, dpi = 400, bg = "white")

  # Compare residual ordination inertia with out-of-bag prediction of genomic axes.
  axis_variance <- RDAforest_ordination$CA$eig /
    sum(RDAforest_ordination$CA$eig)
  axis_variance_table <- data.frame(
    axis = seq_along(axis_variance),
    proportion = axis_variance
  )
  plot_axis_variance <- ggplot(
    axis_variance_table,
    aes(axis, proportion)
  ) +
    geom_segment(
      aes(xend = axis, y = 0, yend = proportion),
      linewidth = .75, colour = idv_purple
    ) +
    geom_vline(
      xintercept = RDAforest_top_pcs + .5,
      linetype = 2, colour = idv_orange, linewidth = .8
    ) +
    theme_ameiva() +
    labs(
      title = "Conditioned genomic ordination",
      x = "Genomic PCoA axis",
      y = "Proportion of residual inertia"
    )

  plot_axis_prediction <- ggplot(
    RDAforest_PC_prediction,
    aes(PC, OOB_R2, fill = OOB_R2)
  ) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = .45) +
    geom_col(width = .72) +
    scale_fill_gradient2(
      low = idv_orange, mid = idv_cream, high = idv_purple,
      midpoint = 0, guide = "none"
    ) +
    theme_ameiva() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(
      title = "Predictability of genomic axes",
      x = NULL, y = "Cross-validated R²"
    )

  png(file.path(figures, "15_RDAforest_axes_predictability.png"),
      width = 2400, height = 1200, res = 220)
  grid::grid.newpage()
  print(plot_axis_variance, vp = grid::viewport(
    x = .25, y = .5, width = .48, height = .94
  ))
  print(plot_axis_prediction, vp = grid::viewport(
    x = .75, y = .5, width = .48, height = .94
  ))
  dev.off()

  plot_mtry_delta <- ggplot(mtry_selection$delta, aes(var, values)) +
    geom_hline(yintercept = 0, colour = idv_orange, linewidth = .8) +
    geom_boxplot(
      outlier.shape = NA, fill = idv_purple_light,
      colour = idv_purple_dark, alpha = .72
    ) +
    coord_flip() + theme_ameiva() +
    labs(
      title = "RDAforest mtry selection",
      subtitle = "Change in importance at higher mtry across 31 jackknife replicates",
      x = NULL, y = "Importance(high mtry) - importance(low mtry)"
    )
  ggsave(file.path(figures, "16_RDAforest_mtry_delta.png"), plot_mtry_delta,
         width = 8, height = 6, dpi = 400, bg = "white")

  plot_mtry_positive <- ggplot(
    mtry_selection$prop.positive,
    aes(var, prop.positive, fill = var %in% selected_RDAforest)
  ) +
    geom_col(width = .72) +
    geom_hline(yintercept = .50, colour = idv_orange_dark, linewidth = .8) +
    coord_flip() + theme_ameiva() +
    scale_fill_manual(values = c(
      "FALSE" = idv_orange_light, "TRUE" = idv_purple
    ),
                      guide = "none") +
    labs(
      title = "Variables retained by the RDAforest mtry criterion",
      subtitle = "Purple variables also exceed 10% of the strongest median importance",
      x = NULL, y = "Proportion of non-negative changes"
    )
  ggsave(file.path(figures, "17_RDAforest_mtry_positive.png"), plot_mtry_positive,
         width = 8, height = 6, dpi = 400, bg = "white")

  plot_importance_uncertainty <- ggplot(
    rdaforest_fit$all.importances,
    aes(variable, importance)
  ) +
    geom_boxplot(
      outlier.shape = NA, fill = idv_purple_light,
      colour = idv_purple_dark, alpha = .70
    ) +
    stat_summary(
      fun = median, geom = "point", shape = 18,
      size = 3, colour = idv_orange
    ) +
    coord_flip() + theme_ameiva() +
    labs(
      title = "RDAforest ordination-jackknife uncertainty",
      subtitle = "25 replicates; 80% of individuals retained per ordination",
      x = NULL, y = "Proportion of total genomic variation explained"
    )
  ggsave(file.path(figures, "18_RDAforest_importance_uncertainty.png"),
         plot_importance_uncertainty, width = 8, height = 6,
         dpi = 400, bg = "white")

  exploratory_variables <- names(sort(
    mtry_selection$importances, decreasing = TRUE
  ))
  exploratory_turnover <- do.call(rbind, lapply(
    exploratory_variables,
    function(variable) {
      curve <- gradientForest::cumimp(
        RDAforest_exploratory, variable
      )
      data.frame(
        environment = variable,
        predictor = curve$x,
        cumulative_importance = curve$y
      )
    }
  ))
  exploratory_turnover$environment <- factor(
    exploratory_turnover$environment,
    levels = exploratory_variables
  )
  plot_exploratory_turnover <- ggplot(
    exploratory_turnover,
    aes(predictor, cumulative_importance)
  ) +
    geom_area(fill = idv_orange_light, alpha = .48) +
    geom_line(colour = idv_purple, linewidth = .95) +
    facet_wrap(~environment, scales = "free_x", ncol = 2) +
    theme_ameiva() +
    labs(
      title = "Exploratory RDAforest turnover",
      subtitle = "Cumulative importance from the initial makeGF model",
      x = "Standardized environmental predictor",
      y = "Cumulative importance"
    )
  ggsave(
    file.path(figures, "19_RDAforest_exploratory_turnover.png"),
    plot_exploratory_turnover, width = 10, height = 7.5,
    dpi = 400, bg = "white"
  )

  jackknife_turnover <- do.call(rbind, lapply(
    selected_RDAforest,
    function(variable) {
      ordered <- order(RDAforest_environment[[variable]])
      data.frame(
        environment = variable,
        predictor = RDAforest_environment[[variable]][ordered],
        cumulative_R2 =
          rdaforest_fit$predictions.turnover[ordered, variable]
      )
    }
  ))
  plot_jackknife_turnover <- ggplot(
    jackknife_turnover,
    aes(predictor, cumulative_R2)
  ) +
    geom_step(colour = idv_purple, linewidth = 1.05) +
    geom_rug(colour = idv_orange, alpha = .28, sides = "b") +
    facet_wrap(~environment, scales = "free_x", nrow = 1) +
    theme_ameiva() +
    labs(
      title = "Jackknifed RDAforest turnover",
      subtitle = "Final predictors retained after mtry selection",
      x = "Standardized environmental predictor",
      y = "Cumulative R²"
    )
  ggsave(
    file.path(figures, "20_RDAforest_jackknife_turnover.png"),
    plot_jackknife_turnover, width = 10, height = 5.5,
    dpi = 400, bg = "white"
  )

  grDevices::pdf(NULL)
  RDAforest_sample_plot <- RDAforest::plot_adaptation(
    rdaforest_fit$predictions.direct,
    RDAforest_environment[, selected_RDAforest, drop = FALSE],
    Variables[, c("long", "lat")],
    main = "Sample-level RDAforest predictions",
    rangeExp = 1.5, scal = 10, jitscale = .05,
    color.scheme = "001", lighten = .65
  )
  dev.off()
  RDAforest_sample_panels <- adaptation_ggplot_panels(
    RDAforest_sample_plot,
    RDAforest_environment[, selected_RDAforest, drop = FALSE]
  )
  save_ggplot_row(
    RDAforest_sample_panels,
    file.path(figures, "21_RDAforest_sample_adaptation.png"),
    width = 2400, height = 1200, res = 220
  )
}


# 17. Describe multilocus genotype clusters among shared candidates ----

mosaic_loci <- candidates_consensus
mosaic_genotypes <- geno_complete[, mosaic_loci, drop = FALSE]
stopifnot(ncol(mosaic_genotypes) >= 2)

mosaic_axes <- min(15, nrow(mosaic_genotypes) - 1, ncol(mosaic_genotypes))
mosaic_pca <- prcomp(
  mosaic_genotypes,
  center = TRUE,
  scale. = FALSE,
  rank. = mosaic_axes
)
mosaic_variance <- mosaic_pca$sdev^2 / sum(mosaic_pca$sdev^2)
mosaic_scores <- sweep(
  mosaic_pca$x[, seq_len(mosaic_axes), drop = FALSE],
  2,
  sqrt(mosaic_variance[seq_len(mosaic_axes)]),
  "*"
)

k_values <- 2:6
kmeans_models <- vector("list", length(k_values))
silhouette_values <- numeric(length(k_values))
bootstrap_ARI <- numeric(length(k_values))

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  observed <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  column_pairs <- sum(choose2(colSums(tab)))
  expected <- row_pairs * column_pairs / choose2(n)
  (observed - expected) /
    ((row_pairs + column_pairs) / 2 - expected)
}

for (i in seq_along(k_values)) {
  set.seed(2027)
  kmeans_models[[i]] <- kmeans(mosaic_scores, centers = k_values[i], nstart = 100)
  silhouette_values[i] <- mean(silhouette(
    kmeans_models[[i]]$cluster,
    dist(mosaic_scores)
  )[, 3])
  bootstrap_ARI[i] <- mean(replicate(200, {
    keep <- sample(
      seq_len(nrow(mosaic_scores)),
      ceiling(.8 * nrow(mosaic_scores))
    )
    subset_fit <- kmeans(
      mosaic_scores[keep, , drop = FALSE],
      centers = k_values[i], nstart = 50
    )
    adjusted_rand_index(
      kmeans_models[[i]]$cluster[keep], subset_fit$cluster
    )
  }))
}

best_mosaic_index <- which.max(silhouette_values)
best_mosaic_k <- k_values[best_mosaic_index]
mosaic_cluster <- factor(kmeans_models[[best_mosaic_index]]$cluster)

mosaic_table <- data.frame(
  InternalID = rownames(geno_complete),
  longitude = Variables$long,
  latitude = Variables$lat,
  PC1 = mosaic_pca$x[, 1],
  PC2 = mosaic_pca$x[, 2],
  genotype_cluster = mosaic_cluster
)

write.csv(mosaic_table, file.path(out, "13_adaptive_mosaic_genotype_clusters.csv"), row.names = FALSE)
write.csv(
  data.frame(
    K = k_values,
    mean_silhouette = silhouette_values,
    bootstrap_adjusted_Rand = bootstrap_ARI,
    selected = k_values == best_mosaic_k
  ),
  file.path(out, "13_adaptive_mosaic_K_diagnostics.csv"),
  row.names = FALSE
)

mosaic_colours <- setNames(
  idv_discrete(best_mosaic_k), levels(mosaic_cluster)
)

plot_mosaic_PCA <- ggplot(mosaic_table, aes(PC1, PC2, colour = genotype_cluster)) +
  geom_point(size = 2.7, alpha = .9) +
  scale_colour_manual(values = mosaic_colours) +
  theme_ameiva() +
  labs(
    title = "Adaptive genotypes form a mosaic within a generalist species",
    colour = "Genotype cluster"
  )

plot_mosaic_map <- ggplot(mosaic_table, aes(longitude, latitude, fill = genotype_cluster)) +
  geom_polygon(data = ameiva_world_map, aes(long, lat, group = group), inherit.aes = FALSE, fill = idv_cream, colour = "grey75", linewidth = .2) +
  geom_point(shape = 21, colour = "white", stroke = .45, size = 3,
             position = position_jitter(width = .12, height = .12, seed = 2027)) +
  scale_fill_manual(values = mosaic_colours) +
  coord_quickmap(xlim = c(-82, -32), ylim = c(-35, 12)) +
  theme_ameiva() +
  labs(
    title = "Adaptive mosaic across South America",
    subtitle = "Point colours are clusters inferred directly from candidate genotypes",
    fill = "Genotype cluster"
  )

ggsave(file.path(figures, "05_adaptive_mosaic_genotype_PCA.png"), plot_mosaic_PCA,
       width = 7, height = 5.5, dpi = 400, bg = "white")
ggsave(file.path(figures, "06_adaptive_mosaic_map.png"), plot_mosaic_map,
       width = 8, height = 7, dpi = 400, bg = "white")

# Compare K = 3 with the K = 2 solution favoured by mean silhouette width.
mosaic_table_K3 <- mosaic_table
mosaic_table_K3$genotype_cluster <- factor(kmeans_models[[which(k_values == 3)]]$cluster)
plot_mosaic_PCA_K3 <- ggplot(
  mosaic_table_K3, aes(PC1, PC2, colour = genotype_cluster)
) +
  geom_point(size = 2.7, alpha = .9) +
  scale_colour_manual(values = idv_discrete(3)) +
  theme_ameiva() +
  labs(
    title = "Candidate-genotype sensitivity at K = 3",
    subtitle = "Displayed because PC2 suggests subdivision; K = 2 has the stronger silhouette",
    colour = "K = 3 cluster"
  )
ggsave(
  file.path(figures, "06b_adaptive_mosaic_genotype_PCA_K3_sensitivity.png"),
  plot_mosaic_PCA_K3, width = 7, height = 5.5, dpi = 400, bg = "white"
)


# 18. Project candidate-panel GF and RDA onto climatic rasters ----

spatial_out <- file.path(out, "spatial_offset")
spatial_temp <- file.path(spatial_out, "terra_tmp")
dir.create(spatial_out, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_temp, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = spatial_temp, progress = 1)

rasters_dir <- "/projects/rstudio-export/Ameiva"
study_extent <- ext(-85, -30, -35, 15)
raster_aggregate_factor <- 5
land_mask <- terra::vect(sf::st_as_sf(maps::map(
  "world", fill = TRUE, plot = FALSE
)))

raster_file_for <- function(bio_name, period) {
  file.path(rasters_dir, paste0("CHELSA_", bio_name, "_", period, "_SA.tif"))
}

read_climate_period <- function(period) {
  raster_files <- vapply(
    names(Environment), raster_file_for, character(1), period = period
  )
  stopifnot(all(file.exists(raster_files)))
  raster_layers <- lapply(names(Environment), function(variable) {
    x <- crop(rast(raster_files[[variable]]), study_extent)
    if (raster_aggregate_factor > 1) {
      x <- terra::aggregate(x, fact = raster_aggregate_factor, fun = mean)
    }
    terra::crs(land_mask) <- terra::crs(x)
    x <- terra::mask(x, land_mask)
    (x - environment_center[variable]) / environment_scale[variable]
  })
  result <- rast(raster_layers)
  names(result) <- names(Environment)
  result
}

if (RUN_SPATIAL_PROJECTION) {

  climate_present <- read_climate_period("1981-2010")
  climate_mid <- read_climate_period("2041-2070")
  climate_late <- read_climate_period("2071-2100")

  # Classify cells by their distance from the sampled climatic ranges.
  # Allow a buffer of 10% of the observed span for present climate and 20% for future climate.
  # Cells beyond either period's allowance are excluded from offset comparisons.
  climate_training_ranges <- vapply(
    Environment,
    range,
    numeric(2)
  )
  rownames(climate_training_ranges) <- c("minimum", "maximum")

  audit_climate_domain <- function(climate_stack, extra_fraction) {
    training_span <-
      climate_training_ranges["maximum", ] -
      climate_training_ranges["minimum", ]
    outside_training <- climate_stack[[1]] * 0
    outside_allowed <- climate_stack[[1]] * 0
    for (variable in names(Environment)) {
      outside_training <- outside_training +
        (climate_stack[[variable]] < climate_training_ranges["minimum", variable] |
         climate_stack[[variable]] > climate_training_ranges["maximum", variable])
      outside_allowed <- outside_allowed +
        (climate_stack[[variable]] <
           climate_training_ranges["minimum", variable] -
           extra_fraction * training_span[variable] |
         climate_stack[[variable]] >
           climate_training_ranges["maximum", variable] +
           extra_fraction * training_span[variable])
    }
    status <- ifel(
      outside_allowed > 0, 3,
      ifel(outside_training > 0, 2, 1)
    )
    names(status) <- "projection_status"
    allowed_mask <- ifel(status < 3, 1, NA)
    names(allowed_mask) <- "allowed_projection"
    list(
      status = status,
      allowed_mask = allowed_mask,
      outside_training = outside_training,
      outside_allowed = outside_allowed
    )
  }

  combine_offset_domain <- function(present_domain, future_domain) {
    status <- ifel(
      present_domain$status == 3 | future_domain$status == 3, 3,
      ifel(
        present_domain$status == 2 | future_domain$status == 2,
        2, 1
      )
    )
    names(status) <- "projection_status"
    allowed_mask <- ifel(status < 3, 1, NA)
    names(allowed_mask) <- "allowed_projection"
    list(status = status, allowed_mask = allowed_mask)
  }

  climate_domain_present <- audit_climate_domain(climate_present, .10)
  climate_domain_mid <- audit_climate_domain(climate_mid, .20)
  climate_domain_late <- audit_climate_domain(climate_late, .20)
  offset_domain_mid <- combine_offset_domain(
    climate_domain_present, climate_domain_mid
  )
  offset_domain_late <- combine_offset_domain(
    climate_domain_present, climate_domain_late
  )

  summarize_projection_domain <- function(domain_status, period) {
    counts <- vapply(1:3, function(code) {
      as.numeric(global(domain_status == code, "sum", na.rm = TRUE)[1, 1])
    }, numeric(1))
    data.frame(
      period = period,
      status = c("inside", "near-range", "excluded"),
      cells = counts,
      proportion = counts / sum(counts)
    )
  }
  offset_domain_summary <- rbind(
    summarize_projection_domain(offset_domain_mid$status, "2041-2070"),
    summarize_projection_domain(offset_domain_late$status, "2071-2100")
  )
  write.csv(
    offset_domain_summary,
    file.path(spatial_out, "19_offset_projection_domain_summary.csv"),
    row.names = FALSE
  )
  writeRaster(
    offset_domain_mid$status,
    file.path(spatial_out, "19_offset_projection_domain_mid.tif"),
    overwrite = TRUE
  )
  writeRaster(
    offset_domain_late$status,
    file.path(spatial_out, "19_offset_projection_domain_late.tif"),
    overwrite = TRUE
  )

  # Compare raster values extracted at sampling coordinates with the original climate table.
  climate_alignment_audit <- do.call(rbind, lapply(
    names(Environment),
    function(variable) {
      raster_raw_units <-
        climate_present[[variable]] * environment_scale[variable] +
        environment_center[variable]
      extracted <- terra::extract(
        raster_raw_units,
        Variables[, c("long", "lat")],
        ID = FALSE
      )[[1]]
      csv_column <- paste0(
        "X.CHELSA_", variable, "_1981.2010_SA"
      )
      observed <- Variables[[csv_column]]
      complete <- is.finite(extracted) & is.finite(observed)
      fit <- lm(extracted[complete] ~ observed[complete])
      data.frame(
        variable = variable,
        n = sum(complete),
        correlation = cor(extracted[complete], observed[complete]),
        intercept = unname(coef(fit)[1]),
        slope = unname(coef(fit)[2]),
        normalized_RMSE =
          sqrt(mean((extracted[complete] - observed[complete])^2)) /
          sd(observed[complete]),
        raster_resolution_degrees = terra::res(climate_present)[1],
        geometry_identical_across_periods =
          compareGeom(
            climate_present[[variable]], climate_mid[[variable]],
            climate_late[[variable]], stopOnError = FALSE
          )
      )
    }
  ))
  write.csv(
    climate_alignment_audit,
    file.path(spatial_out, "19_climate_input_alignment_audit.csv"),
    row.names = FALSE
  )

  # Transform each predictor through its GF cumulative-importance curve.
  # Fit the display PCA to transformed values at sampled individuals.
  GF_transformed_training <- as.data.frame(predict(
    gradient_forest, newdata = gradient_forest$X, extrap = FALSE
  ))
  GF_turnover_pca <- prcomp(
    GF_transformed_training[, names(Environment), drop = FALSE],
    center = TRUE,
    scale. = TRUE
  )

  transform_GF_stack <- function(climate_stack) {
    transformed_layers <- lapply(names(Environment), function(variable) {
      curve <- gradientForest::cumimp(gradient_forest, variable)
      app(climate_stack[[variable]], function(values) {
        approx(curve$x, curve$y, xout = values, rule = 2)$y
      })
    })
    result <- rast(transformed_layers)
    names(result) <- names(Environment)
    result
  }

  GF_present_transformed <- transform_GF_stack(climate_present)
  GF_mid_transformed <- transform_GF_stack(climate_mid)
  GF_late_transformed <- transform_GF_stack(climate_late)

  project_GF_PCA <- function(transformed_stack) {
    standardized <- transformed_stack
    for (variable in names(Environment)) {
      standardized[[variable]] <- (
        transformed_stack[[variable]] - GF_turnover_pca$center[variable]
      ) / GF_turnover_pca$scale[variable]
    }
    axes <- lapply(1:3, function(axis) {
      sum(standardized * GF_turnover_pca$rotation[names(Environment), axis])
    })
    result <- rast(axes)
    names(result) <- paste0("GF_PC", 1:3)
    result
  }

  GF_present_PCA <- project_GF_PCA(GF_present_transformed)
  GF_mid_PCA <- project_GF_PCA(GF_mid_transformed)
  GF_late_PCA <- project_GF_PCA(GF_late_transformed)

  # Use present-day empirical percentiles as a shared colour scale across periods.
  set.seed(2027)
  GF_ECDF_reference <- lapply(1:3, function(axis) {
    reference <- spatSample(
      GF_present_PCA[[axis]], size = 200000,
      method = "random", na.rm = TRUE, as.df = TRUE
    )[[1]]
    sort(reference[is.finite(reference)])
  })
  make_GF_RGB <- function(pca_stack) {
    rgb_layers <- lapply(1:3, function(axis) {
      reference <- GF_ECDF_reference[[axis]]
      probabilities <- seq_along(reference) / length(reference)
      app(pca_stack[[axis]], function(values) {
        approx(
          reference, probabilities, xout = values,
          rule = 2, ties = "ordered"
        )$y
      })
    })
    result <- rast(rgb_layers)
    names(result) <- c("red", "green", "blue")
    result
  }

  GF_present_RGB <- make_GF_RGB(GF_present_PCA)
  GF_mid_RGB <- make_GF_RGB(GF_mid_PCA)
  GF_late_RGB <- make_GF_RGB(GF_late_PCA)

  writeRaster(GF_present_PCA, file.path(spatial_out, "18_GF_present_PCA.tif"), overwrite = TRUE)
  writeRaster(GF_present_RGB, file.path(spatial_out, "18_GF_present_RGB.tif"), overwrite = TRUE)
  writeRaster(GF_mid_RGB, file.path(spatial_out, "18_GF_mid_century_RGB.tif"), overwrite = TRUE)
  writeRaster(GF_late_RGB, file.path(spatial_out, "18_GF_late_century_RGB.tif"), overwrite = TRUE)
  saveRDS(GF_turnover_pca, file.path(spatial_out, "18_GF_turnover_PCA.rds"))

  # Fit a climate-only RDA to the strict landscape candidate subset.
  # The partial RDA was used for detection; this refit provides spatial projection coefficients.
  stopifnot(length(candidates_RDA_material_strict) >= 2)
  RDA_projection_model <- rda(
    geno_complete[, candidates_RDA_material_strict, drop = FALSE] ~ .,
    data = Environment
  )
  RDA_projection_axes <- min(3, RDA_projection_model$CCA$rank)
  RDA_environment_loadings <- RDA_projection_model$CCA$biplot[
    names(Environment), seq_len(RDA_projection_axes), drop = FALSE
  ]
  RDA_axis_weights <- RDA_projection_model$CCA$eig[
    seq_len(RDA_projection_axes)
  ] / sum(RDA_projection_model$CCA$eig[seq_len(RDA_projection_axes)])

  project_RDA <- function(climate_stack) {
    axes <- lapply(seq_len(RDA_projection_axes), function(axis) {
      sum(climate_stack * RDA_environment_loadings[, axis])
    })
    result <- rast(axes)
    names(result) <- paste0("RDA", seq_len(RDA_projection_axes))
    result
  }

  RDA_present <- project_RDA(climate_present)
  RDA_mid <- project_RDA(climate_mid)
  RDA_late <- project_RDA(climate_late)
  writeRaster(RDA_present, file.path(spatial_out, "18_RDA_present_axes.tif"), overwrite = TRUE)
  saveRDS(RDA_projection_model, file.path(spatial_out, "18_RDA_projection_model.rds"))


# 19. Calculate candidate-panel GF and RDA genomic offset ----

  GF_offset_mid <-
    sqrt(sum((GF_mid_transformed - GF_present_transformed)^2)) *
    offset_domain_mid$allowed_mask
  GF_offset_late <-
    sqrt(sum((GF_late_transformed - GF_present_transformed)^2)) *
    offset_domain_late$allowed_mask
  names(GF_offset_mid) <- "GF_offset_mid_century"
  names(GF_offset_late) <- "GF_offset_late_century"

  RDA_offset <- function(future_axes) {
    differences <- future_axes - RDA_present
    weighted <- lapply(seq_len(RDA_projection_axes), function(axis) {
      RDA_axis_weights[axis] * differences[[axis]]^2
    })
    sqrt(sum(rast(weighted)))
  }
  RDA_offset_mid <- RDA_offset(RDA_mid) * offset_domain_mid$allowed_mask
  RDA_offset_late <- RDA_offset(RDA_late) * offset_domain_late$allowed_mask
  names(RDA_offset_mid) <- "RDA_offset_mid_century"
  names(RDA_offset_late) <- "RDA_offset_late_century"

  writeRaster(GF_offset_mid, file.path(spatial_out, "19_GF_offset_mid_century.tif"), overwrite = TRUE)
  writeRaster(GF_offset_late, file.path(spatial_out, "19_GF_offset_late_century.tif"), overwrite = TRUE)
  writeRaster(RDA_offset_mid, file.path(spatial_out, "19_RDA_offset_mid_century.tif"), overwrite = TRUE)
  writeRaster(RDA_offset_late, file.path(spatial_out, "19_RDA_offset_late_century.tif"), overwrite = TRUE)

  normalize_raster <- function(x) {
    limits <- as.numeric(global(x, range, na.rm = TRUE)[1,])
    (x - limits[1]) / diff(limits)
  }
  disagreement_mid <- abs(normalize_raster(RDA_offset_mid) - normalize_raster(GF_offset_mid))
  disagreement_late <- abs(normalize_raster(RDA_offset_late) - normalize_raster(GF_offset_late))
  names(disagreement_mid) <- "RDA_GF_disagreement_mid"
  names(disagreement_late) <- "RDA_GF_disagreement_late"
  writeRaster(disagreement_mid, file.path(spatial_out, "19_RDA_GF_disagreement_mid.tif"), overwrite = TRUE)
  writeRaster(disagreement_late, file.path(spatial_out, "19_RDA_GF_disagreement_late.tif"), overwrite = TRUE)

  set.seed(2030)
  offset_sample_mid <- spatSample(
    c(RDA_offset_mid, GF_offset_mid), 50000,
    method = "random", na.rm = TRUE, as.df = TRUE
  )
  set.seed(2031)
  offset_sample_late <- spatSample(
    c(RDA_offset_late, GF_offset_late), 50000,
    method = "random", na.rm = TRUE, as.df = TRUE
  )
  offset_summary <- data.frame(
    period = c("2041-2070", "2071-2100"),
    RDA_mean = c(mean(offset_sample_mid[[1]]), mean(offset_sample_late[[1]])),
    GF_mean = c(mean(offset_sample_mid[[2]]), mean(offset_sample_late[[2]])),
    RDA_GF_spearman = c(
      cor(offset_sample_mid[[1]], offset_sample_mid[[2]], method = "spearman"),
      cor(offset_sample_late[[1]], offset_sample_late[[2]], method = "spearman")
    )
  )
  write.csv(offset_summary, file.path(spatial_out, "19_offset_summary.csv"), row.names = FALSE)


# 20. Classify projected landscapes and visualise genomic offset ----

  # Aggregate offset maps to approximately 0.5 degrees for display; numerical
  # analyses retain the 2.5-arc-minute grid.
  offset_display_factor <- max(
    1,
    round(.5 / terra::res(GF_offset_mid)[1])
  )
  GF_offset_mid_display <- aggregate(
    GF_offset_mid, fact = offset_display_factor,
    fun = mean, na.rm = TRUE
  )
  GF_offset_late_display <- aggregate(
    GF_offset_late, fact = offset_display_factor,
    fun = mean, na.rm = TRUE
  )
  RDA_offset_mid_display <- aggregate(
    RDA_offset_mid, fact = offset_display_factor,
    fun = mean, na.rm = TRUE
  )
  RDA_offset_late_display <- aggregate(
    RDA_offset_late, fact = offset_display_factor,
    fun = mean, na.rm = TRUE
  )
  disagreement_mid_display <- aggregate(
    disagreement_mid, fact = offset_display_factor,
    fun = mean, na.rm = TRUE
  )
  disagreement_late_display <- aggregate(
    disagreement_late, fact = offset_display_factor,
    fun = mean, na.rm = TRUE
  )
  domain_mid_display <- aggregate(
    offset_domain_mid$status, fact = offset_display_factor,
    fun = modal, na.rm = TRUE
  )
  domain_late_display <- aggregate(
    offset_domain_late$status, fact = offset_display_factor,
    fun = modal, na.rm = TRUE
  )
  excluded_mid_display <- ifel(domain_mid_display == 3, 1, NA)
  excluded_late_display <- ifel(domain_late_display == 3, 1, NA)

  classify_landscape <- function(score_stack, seed = 2027) {
    set.seed(seed)
    landscape <- spatSample(
      score_stack, size = 50000, method = "random",
      na.rm = TRUE, as.df = TRUE
    )
    landscape <- landscape[complete.cases(landscape),,drop = FALSE]
    landscape_center <- vapply(landscape, mean, numeric(1))
    landscape_scale <- vapply(landscape, sd, numeric(1))
    landscape_z <- scale(landscape, center = landscape_center, scale = landscape_scale)
    set.seed(seed)
    diagnostic_rows <- sample(seq_len(nrow(landscape_z)), min(4000, nrow(landscape_z)))
    k_values <- 2:6
    silhouette_mean <- vapply(k_values, function(k) {
      fit_k <- kmeans(landscape_z[diagnostic_rows,,drop = FALSE], k, nstart = 30)
      mean(silhouette(fit_k$cluster, dist(landscape_z[diagnostic_rows,,drop = FALSE]))[,3])
    }, numeric(1))
    selected_k <- k_values[which.max(silhouette_mean)]
    set.seed(seed)
    fit <- kmeans(landscape_z, selected_k, nstart = 50)

    standardized_stack <- score_stack
    for (axis in seq_len(nlyr(score_stack))) {
      standardized_stack[[axis]] <- (
        score_stack[[axis]] - landscape_center[axis]
      ) / landscape_scale[axis]
    }
    distances <- rast(lapply(seq_len(selected_k), function(k) {
      sum((standardized_stack - fit$centers[k,])^2)
    }))
    regions <- distances[[1]] * 0 + 1
    minimum_distance <- distances[[1]]
    if (selected_k > 1) for (k in 2:selected_k) {
      closer <- distances[[k]] < minimum_distance
      regions <- ifel(closer, k, regions)
      minimum_distance <- ifel(closer, distances[[k]], minimum_distance)
    }
    list(
      regions = regions,
      selected_k = selected_k,
      diagnostics = data.frame(K = k_values, silhouette = silhouette_mean),
      fit = fit
    )
  }

  GF_landscape <- classify_landscape(GF_present_PCA, seed = 2027)
  RDA_landscape <- classify_landscape(RDA_present, seed = 2028)
  names(GF_landscape$regions) <- "GF_landscape_cluster"
  names(RDA_landscape$regions) <- "RDA_landscape_cluster"
  writeRaster(GF_landscape$regions, file.path(spatial_out, "20_GF_present_kmeans.tif"), overwrite = TRUE)
  writeRaster(RDA_landscape$regions, file.path(spatial_out, "20_RDA_present_kmeans.tif"), overwrite = TRUE)
  write.csv(
    rbind(
      transform(GF_landscape$diagnostics, method = "Gradient Forest"),
      transform(RDA_landscape$diagnostics, method = "RDA")
    ),
    file.path(spatial_out, "20_landscape_kmeans_diagnostics.csv"),
    row.names = FALSE
  )

  map_display_factor <- max(
    1, round(.25 / terra::res(GF_present_RGB)[1])
  )
  GF_present_RGB_display <- aggregate(
    GF_present_RGB, fact = map_display_factor,
    fun = mean, na.rm = TRUE
  )
  GF_mid_RGB_display <- aggregate(
    GF_mid_RGB, fact = map_display_factor,
    fun = mean, na.rm = TRUE
  )
  GF_late_RGB_display <- aggregate(
    GF_late_RGB, fact = map_display_factor,
    fun = mean, na.rm = TRUE
  )

  GF_present_colour_data <- as.data.frame(
    GF_present_RGB_display, xy = TRUE, na.rm = TRUE
  )
  GF_present_colour_data$brand_colour <- idv_false_colour(
    GF_present_colour_data$red,
    GF_present_colour_data$green,
    GF_present_colour_data$blue
  )
  sample_map_data <- data.frame(
    longitude = Variables$long,
    latitude = Variables$lat,
    genotype_cluster = mosaic_cluster
  )
  plot_GF_present <- ggplot() +
    geom_polygon(data = ameiva_world_map, aes(long, lat, group = group), inherit.aes = FALSE, fill = idv_cream, colour = "grey78", linewidth = .2) +
    geom_raster(
      data = GF_present_colour_data,
      aes(x, y, fill = brand_colour)
    ) +
    scale_fill_identity() +
    geom_point(
      data = sample_map_data,
      aes(longitude, latitude, colour = genotype_cluster),
      size = 2.2, alpha = .9
    ) +
    scale_colour_manual(values = mosaic_colours) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Present adaptive landscape — Gradient Forest",
      subtitle = "Three ECDF-transformed turnover axes in purple, orange and teal",
      x = "Longitude", y = "Latitude",
      colour = "Genotype cluster"
    )
  ggsave(
    file.path(figures, "07_GF_present_RGB.png"),
    plot_GF_present, width = 8, height = 6,
    dpi = 400, bg = "white"
  )

  GF_future_colour_data <- rbind(
    transform(
      as.data.frame(GF_mid_RGB_display, xy = TRUE, na.rm = TRUE),
      period = "2041-2070"
    ),
    transform(
      as.data.frame(GF_late_RGB_display, xy = TRUE, na.rm = TRUE),
      period = "2071-2100"
    )
  )
  GF_future_colour_data$brand_colour <- idv_false_colour(
    GF_future_colour_data$red,
    GF_future_colour_data$green,
    GF_future_colour_data$blue
  )
  plot_GF_future <- ggplot(
    GF_future_colour_data,
    aes(x, y, fill = brand_colour)
  ) +
    geom_raster() +
    scale_fill_identity() +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap() +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Future Gradient Forest turnover",
      subtitle = "Same present-defined ECDF colour scale in both periods",
      x = "Longitude", y = "Latitude"
    )
  ggsave(
    file.path(figures, "08_GF_future_RGB.png"),
    plot_GF_future, width = 10, height = 5.5,
    dpi = 400, bg = "white"
  )

  excluded_plot_data <- rbind(
    transform(
      raster_to_ggplot(excluded_mid_display, "excluded"),
      period = "2041-2070"
    ),
    transform(
      raster_to_ggplot(excluded_late_display, "excluded"),
      period = "2071-2100"
    )
  )
  excluded_plot_data <- excluded_plot_data[
    is.finite(excluded_plot_data$excluded), , drop = FALSE
  ]

  GF_offset_plot_data <- rbind(
    raster_to_ggplot(
      GF_offset_mid_display, "offset", "2041-2070"
    ),
    raster_to_ggplot(
      GF_offset_late_display, "offset", "2071-2100"
    )
  )
  GF_offset_plot_data <- GF_offset_plot_data[
    is.finite(GF_offset_plot_data$offset), , drop = FALSE
  ]
  plot_GF_offset <- ggplot() +
    geom_raster(
      data = excluded_plot_data,
      aes(x, y), fill = idv_grey
    ) +
    geom_raster(
      data = GF_offset_plot_data,
      aes(x, y, fill = offset)
    ) +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap() +
    scale_fill_gradientn(colours = idv_sequential) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Gradient Forest genomic offset",
      subtitle = "Grey cells exceed the allowed environmental domain",
      x = "Longitude", y = "Latitude", fill = "GF offset"
    )
  ggsave(
    file.path(figures, "09_GF_genomic_offset.png"),
    plot_GF_offset, width = 10, height = 5.5,
    dpi = 400, bg = "white"
  )

  RDA_present_display <- aggregate(
    RDA_present[[1]], fact = map_display_factor,
    fun = mean, na.rm = TRUE
  )
  RDA_present_plot_data <- raster_to_ggplot(
    RDA_present_display, "RDA1"
  )
  RDA_present_plot_data <- RDA_present_plot_data[
    is.finite(RDA_present_plot_data$RDA1), , drop = FALSE
  ]
  plot_RDA_present <- ggplot() +
    geom_polygon(data = ameiva_world_map, aes(long, lat, group = group), inherit.aes = FALSE, fill = idv_cream, colour = "grey78", linewidth = .2) +
    geom_raster(
      data = RDA_present_plot_data,
      aes(x, y, fill = RDA1)
    ) +
    geom_point(
      data = sample_map_data,
      aes(longitude, latitude, colour = genotype_cluster),
      size = 2.2, alpha = .9
    ) +
    scale_fill_gradientn(colours = idv_diverging) +
    scale_colour_manual(values = mosaic_colours) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Present adaptive landscape — RDA1",
      x = "Longitude", y = "Latitude",
      fill = "RDA1", colour = "Genotype cluster"
    )
  ggsave(
    file.path(figures, "10_RDA_present_landscape.png"),
    plot_RDA_present, width = 8, height = 6,
    dpi = 400, bg = "white"
  )

  RDA_offset_plot_data <- rbind(
    raster_to_ggplot(
      RDA_offset_mid_display, "offset", "2041-2070"
    ),
    raster_to_ggplot(
      RDA_offset_late_display, "offset", "2071-2100"
    )
  )
  RDA_offset_plot_data <- RDA_offset_plot_data[
    is.finite(RDA_offset_plot_data$offset), , drop = FALSE
  ]
  plot_RDA_offset <- ggplot() +
    geom_raster(
      data = excluded_plot_data,
      aes(x, y), fill = idv_grey
    ) +
    geom_raster(
      data = RDA_offset_plot_data,
      aes(x, y, fill = offset)
    ) +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap() +
    scale_fill_gradientn(colours = idv_sequential) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Eigenvalue-weighted RDA genomic offset",
      subtitle = "Grey cells exceed the allowed environmental domain",
      x = "Longitude", y = "Latitude", fill = "RDA offset"
    )
  ggsave(
    file.path(figures, "11_RDA_genomic_offset.png"),
    plot_RDA_offset, width = 10, height = 5.5,
    dpi = 400, bg = "white"
  )

  disagreement_plot_data <- rbind(
    raster_to_ggplot(
      disagreement_mid_display, "disagreement", "2041-2070"
    ),
    raster_to_ggplot(
      disagreement_late_display, "disagreement", "2071-2100"
    )
  )
  disagreement_plot_data <- disagreement_plot_data[
    is.finite(disagreement_plot_data$disagreement), , drop = FALSE
  ]
  plot_disagreement <- ggplot() +
    geom_raster(
      data = excluded_plot_data,
      aes(x, y), fill = idv_grey
    ) +
    geom_raster(
      data = disagreement_plot_data,
      aes(x, y, fill = disagreement)
    ) +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap() +
    scale_fill_gradientn(colours = idv_sequential) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "RDA–Gradient Forest offset disagreement",
      subtitle = "Absolute difference after independent 0–1 normalization",
      x = "Longitude", y = "Latitude", fill = "Disagreement"
    )
  ggsave(
    file.path(figures, "12_RDA_GF_offset_disagreement.png"),
    plot_disagreement, width = 10, height = 5.5,
    dpi = 400, bg = "white"
  )

  RDA_regions_display <- aggregate(
    RDA_landscape$regions, fact = map_display_factor,
    fun = modal, na.rm = TRUE
  )
  GF_regions_display <- aggregate(
    GF_landscape$regions, fact = map_display_factor,
    fun = modal, na.rm = TRUE
  )
  landscape_region_data <- rbind(
    transform(
      raster_to_ggplot(RDA_regions_display, "cluster"),
      method = paste0("RDA (K = ", RDA_landscape$selected_k, ")")
    ),
    transform(
      raster_to_ggplot(GF_regions_display, "cluster"),
      method = paste0("Gradient Forest (K = ",
                      GF_landscape$selected_k, ")")
    )
  )
  landscape_region_data <- landscape_region_data[
    is.finite(landscape_region_data$cluster), , drop = FALSE
  ]
  maximum_landscape_k <- max(
    RDA_landscape$selected_k, GF_landscape$selected_k
  )
  landscape_region_data$cluster <- factor(
    landscape_region_data$cluster,
    levels = seq_len(maximum_landscape_k)
  )
  plot_landscape_regions <- ggplot(
    landscape_region_data,
    aes(x, y, fill = cluster)
  ) +
    geom_raster() +
    facet_wrap(~method, nrow = 1) +
    coord_quickmap() +
    scale_fill_manual(values = idv_discrete(maximum_landscape_k)) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Adaptive landscape regions",
      subtitle = "Independent k-means summaries of continuous surfaces",
      x = "Longitude", y = "Latitude", fill = "Cluster"
    )
  ggsave(
    file.path(figures, "13_landscape_kmeans_RDA_GF.png"),
    plot_landscape_regions, width = 10, height = 5.5,
    dpi = 400, bg = "white"
  )

  # Display the combined present and future climatic domain.
  domain_plot_data <- rbind(
    transform(
      as.data.frame(domain_mid_display, xy = TRUE, na.rm = TRUE),
      period = "2041-2070"
    ),
    transform(
      as.data.frame(domain_late_display, xy = TRUE, na.rm = TRUE),
      period = "2071-2100"
    )
  )
  domain_plot_data$projection_status <- factor(
    domain_plot_data$projection_status,
    levels = 1:3,
    labels = c("inside", "near-range", "excluded")
  )
  plot_offset_domain <- ggplot(
    domain_plot_data,
    aes(x, y, fill = projection_status)
  ) +
    geom_raster() +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap() +
    scale_fill_manual(values = c(
      inside = idv_purple,
      `near-range` = idv_orange,
      excluded = idv_grey
    )) +
    theme_ameiva() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Projection-domain audit for the four-variable GF/RDA offsets",
      subtitle = paste0(
        "Purple: observed range; orange: within allowance; grey: excluded; ",
        "display grid ≈ ", round(.5, 1), "°"
      ),
      fill = "Projection status"
    )
  ggsave(
    file.path(figures, "29_GF_RDA_offset_domain_audit.png"),
    plot_offset_domain, width = 11, height = 5.5,
    dpi = 400, bg = "white"
  )

  # Inspect future-minus-present climatic changes separately from genomic predictions.
  climate_delta_mid_display <- aggregate(
    climate_mid - climate_present,
    fact = offset_display_factor, fun = mean, na.rm = TRUE
  )
  climate_delta_late_display <- aggregate(
    climate_late - climate_present,
    fact = offset_display_factor, fun = mean, na.rm = TRUE
  )
  climate_delta_plot_data <- do.call(rbind, lapply(
    names(Environment),
    function(variable) {
      paired_delta <- c(
        climate_delta_mid_display[[variable]],
        climate_delta_late_display[[variable]]
      )
      delta_limits <- as.numeric(unlist(global(
        paired_delta, quantile,
        probs = c(.01, .99), na.rm = TRUE
      )))
      symmetric_limit <- max(abs(delta_limits), na.rm = TRUE)
      if (!is.finite(symmetric_limit) || symmetric_limit == 0) {
        symmetric_limit <- 1
      }
      variable_data <- rbind(
        raster_to_ggplot(
          climate_delta_mid_display[[variable]],
          "climate_change", "2041-2070"
        ),
        raster_to_ggplot(
          climate_delta_late_display[[variable]],
          "climate_change", "2071-2100"
        )
      )
      variable_data$environment <- variable
      variable_data$relative_change <- pmax(
        -1, pmin(1, variable_data$climate_change / symmetric_limit)
      )
      variable_data
    }
  ))
  climate_delta_plot_data <- climate_delta_plot_data[
    is.finite(climate_delta_plot_data$relative_change), , drop = FALSE
  ]
  plot_climate_delta_audit <- ggplot(
    climate_delta_plot_data,
    aes(x, y, fill = relative_change)
  ) +
    geom_raster() +
    facet_grid(period ~ environment) +
    coord_quickmap() +
    scale_fill_gradientn(
      colours = idv_diverging,
      limits = c(-1, 1),
      breaks = c(-1, 0, 1),
      labels = c("decrease", "stable", "increase")
    ) +
    theme_ameiva(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Climate-input artifact audit",
      subtitle = paste0(
        "Future minus present; each predictor is symmetrically scaled by its ",
        "paired 1st–99th percentile"
      ),
      x = "Longitude", y = "Latitude",
      fill = "Within-predictor\nrelative change"
    )
  ggsave(
    file.path(figures, "30_climate_change_input_artifact_audit.png"),
    plot_climate_delta_audit, width = 14, height = 7.5,
    dpi = 400, bg = "white"
  )
}


# 21. Project RDAforest and compare genomic offsets ----

# Calculate offset from direct ordination-jackknife predictions using gen_offset_oj(),
# scaled by the 90th percentile of present-day genomic disparity.

RDAforest_spatial_out <- file.path(
  spatial_out,
  if (RDAforest_top_pcs == RDAforest_elbow_pcs) {
    "rdaforest_original"
  } else {
    paste0("rdaforest_pcs", RDAforest_top_pcs, "_geo")
  }
)
dir.create(RDAforest_spatial_out, recursive = TRUE, showWarnings = FALSE)

if (RUN_RDAFOREST && RUN_RDAFOREST_SPATIAL && RUN_SPATIAL_PROJECTION) {

  climate_present_RDAforest <- climate_present[[selected_RDAforest]]
  climate_mid_RDAforest <- climate_mid[[selected_RDAforest]]
  climate_late_RDAforest <- climate_late[[selected_RDAforest]]
  names(climate_present_RDAforest) <- paste0("present_", selected_RDAforest)
  names(climate_mid_RDAforest) <- paste0("mid_", selected_RDAforest)
  names(climate_late_RDAforest) <- paste0("late_", selected_RDAforest)

  # Construct a regular two-dimensional grid and retain cells with complete climate
  # across all three periods. Predict all periods within the same jackknife fit.
  RDAforest_grid_extent <- terra::ext(climate_present_RDAforest)
  RDAforest_grid_aspect <-
    (terra::xmax(RDAforest_grid_extent) - terra::xmin(RDAforest_grid_extent)) /
    (terra::ymax(RDAforest_grid_extent) - terra::ymin(RDAforest_grid_extent))
  RDAforest_grid_nx <- round(sqrt(10000 * RDAforest_grid_aspect))
  RDAforest_grid_ny <- ceiling(10000 / RDAforest_grid_nx)
  RDAforest_grid_dx <-
    (terra::xmax(RDAforest_grid_extent) - terra::xmin(RDAforest_grid_extent)) /
    RDAforest_grid_nx
  RDAforest_grid_dy <-
    (terra::ymax(RDAforest_grid_extent) - terra::ymin(RDAforest_grid_extent)) /
    RDAforest_grid_ny
  RDAforest_spatial_xy_all <- expand.grid(
    x = seq(
      terra::xmin(RDAforest_grid_extent) + RDAforest_grid_dx / 2,
      terra::xmax(RDAforest_grid_extent) - RDAforest_grid_dx / 2,
      length.out = RDAforest_grid_nx
    ),
    y = seq(
      terra::ymin(RDAforest_grid_extent) + RDAforest_grid_dy / 2,
      terra::ymax(RDAforest_grid_extent) - RDAforest_grid_dy / 2,
      length.out = RDAforest_grid_ny
    )
  )
  RDAforest_spatial_values <- terra::extract(
    c(
      climate_present_RDAforest,
      climate_mid_RDAforest,
      climate_late_RDAforest
    ),
    RDAforest_spatial_xy_all,
    ID = FALSE
  )
  RDAforest_spatial_grid <- cbind(
    RDAforest_spatial_xy_all,
    RDAforest_spatial_values
  )
  RDAforest_spatial_grid <- RDAforest_spatial_grid[
    complete.cases(RDAforest_spatial_grid), , drop = FALSE
  ]

  RDAforest_spatial_xy <- RDAforest_spatial_grid[, c("x", "y"), drop = FALSE]
  RDAforest_present_X <- RDAforest_spatial_grid[
    , paste0("present_", selected_RDAforest), drop = FALSE
  ]
  RDAforest_mid_X <- RDAforest_spatial_grid[
    , paste0("mid_", selected_RDAforest), drop = FALSE
  ]
  RDAforest_late_X <- RDAforest_spatial_grid[
    , paste0("late_", selected_RDAforest), drop = FALSE
  ]
  names(RDAforest_present_X) <- selected_RDAforest
  names(RDAforest_mid_X) <- selected_RDAforest
  names(RDAforest_late_X) <- selected_RDAforest

  RDAforest_training_ranges <- apply(
    RDAforest_environment[, selected_RDAforest, drop = FALSE],
    2, range
  )
  audit_RDAforest_domain <- function(new_environment, extra_fraction, period) {
    training_spans <-
      RDAforest_training_ranges[2, ] - RDAforest_training_ranges[1, ]
    outside_training <- matrix(
      FALSE, nrow(new_environment), length(selected_RDAforest),
      dimnames = list(NULL, selected_RDAforest)
    )
    outside_allowed <- outside_training
    for (variable in selected_RDAforest) {
      outside_training[, variable] <-
        new_environment[[variable]] < RDAforest_training_ranges[1, variable] |
        new_environment[[variable]] > RDAforest_training_ranges[2, variable]
      outside_allowed[, variable] <-
        new_environment[[variable]] <
          RDAforest_training_ranges[1, variable] -
          extra_fraction * training_spans[variable] |
        new_environment[[variable]] >
          RDAforest_training_ranges[2, variable] +
          extra_fraction * training_spans[variable]
    }
    status <- ifelse(
      rowSums(outside_allowed) > 0, "excluded",
      ifelse(rowSums(outside_training) > 0, "clamped", "inside")
    )
    data.frame(
      x = RDAforest_spatial_xy$x,
      y = RDAforest_spatial_xy$y,
      period = period,
      status = factor(status, levels = c("inside", "clamped", "excluded")),
      variables_outside_training = rowSums(outside_training),
      variables_outside_allowed_domain = rowSums(outside_allowed)
    )
  }

  RDAforest_domain_audit <- rbind(
    audit_RDAforest_domain(RDAforest_present_X, .10, "Present"),
    audit_RDAforest_domain(RDAforest_mid_X, .20, "2041-2070"),
    audit_RDAforest_domain(RDAforest_late_X, .20, "2071-2100")
  )
  write.csv(
    RDAforest_domain_audit,
    file.path(RDAforest_spatial_out, "21_RDAforest_projection_domain_pixels.csv"),
    row.names = FALSE
  )
  RDAforest_domain_summary <- as.data.frame(table(
    period = RDAforest_domain_audit$period,
    status = RDAforest_domain_audit$status
  ))
  RDAforest_domain_summary$proportion <- ave(
    RDAforest_domain_summary$Freq,
    RDAforest_domain_summary$period,
    FUN = function(x) x / sum(x)
  )
  write.csv(
    RDAforest_domain_summary,
    file.path(RDAforest_spatial_out, "21_RDAforest_projection_domain_summary.csv"),
    row.names = FALSE
  )

  RDAforest_spatial_signature <- data.frame(
    package = as.character(packageVersion("RDAforest")),
    grid_method = "two_dimensional_regular_v1",
    selected_variables = paste(selected_RDAforest, collapse = "+"),
    grid_rows = nrow(RDAforest_spatial_grid),
    combined_prediction_rows = 3 * nrow(RDAforest_spatial_grid),
    top_pcs = RDAforest_top_pcs,
    jackknife_replicates = 25,
    trees = 1500
  )
  RDAforest_signature_file <- file.path(
    RDAforest_spatial_out, "21_RDAforest_spatial_signature.csv"
  )
  RDAforest_combined_file <- file.path(
    RDAforest_spatial_out, "21_RDAforest_combined_ordinationJackknife.rds"
  )
  RDAforest_grid_file <- file.path(
    RDAforest_spatial_out, "21_RDAforest_spatial_grid.rds"
  )

  reuse_RDAforest_spatial <- USE_SAVED_RDS && all(file.exists(c(
    RDAforest_signature_file,
    RDAforest_combined_file,
    RDAforest_grid_file
  ))) && isTRUE(all.equal(
    read.csv(RDAforest_signature_file, check.names = FALSE),
    RDAforest_spatial_signature,
    check.attributes = FALSE
  ))

  if (reuse_RDAforest_spatial) {
    RDAforest_spatial_grid <- readRDS(RDAforest_grid_file)
    RDAforest_spatial_xy <- RDAforest_spatial_grid$xy
    RDAforest_present_X <- RDAforest_spatial_grid$present
    RDAforest_mid_X <- RDAforest_spatial_grid$mid
    RDAforest_late_X <- RDAforest_spatial_grid$late
    RDAforest_combined_spatial <- readRDS(RDAforest_combined_file)
  } else {
    prepare_RDAforest_newX <- function(new_environment, extra_fraction) {
      training_ranges <- apply(
        RDAforest_environment[, selected_RDAforest, drop = FALSE],
        2, range
      )
      training_spans <- training_ranges[2, ] - training_ranges[1, ]
      result <- new_environment
      for (variable in selected_RDAforest) {
        lower_extended <- training_ranges[1, variable] -
          extra_fraction * training_spans[variable]
        upper_extended <- training_ranges[2, variable] +
          extra_fraction * training_spans[variable]
        result[[variable]][result[[variable]] < lower_extended] <- NA
        result[[variable]][result[[variable]] > upper_extended] <- NA
        result[[variable]][result[[variable]] < training_ranges[1, variable]] <-
          training_ranges[1, variable]
        result[[variable]][result[[variable]] > training_ranges[2, variable]] <-
          training_ranges[2, variable]
      }
      result
    }

    RDAforest_newX_combined <- rbind(
      prepare_RDAforest_newX(RDAforest_present_X, extra_fraction = .10),
      prepare_RDAforest_newX(RDAforest_mid_X, extra_fraction = .20),
      prepare_RDAforest_newX(RDAforest_late_X, extra_fraction = .20)
    )

    set.seed(2031)
    RDAforest_combined_spatial <- RDAforest::ordinationJackknife(
      Y = genomic_distance,
      X = RDAforest_environment[, selected_RDAforest, drop = FALSE],
      newX = RDAforest_newX_combined,
      covariates = RDAforest_covariates,
      nreps = 25,
      oob = .20,
      top.pcs = RDAforest_top_pcs,
      ntrees = 1500,
      extra = .20
    )

    saveRDS(
      list(
        xy = RDAforest_spatial_xy,
        present = RDAforest_present_X,
        mid = RDAforest_mid_X,
        late = RDAforest_late_X
      ),
      RDAforest_grid_file
    )
    saveRDS(RDAforest_combined_spatial, RDAforest_combined_file)
    write.csv(
      RDAforest_spatial_signature,
      RDAforest_signature_file,
      row.names = FALSE
    )
  }

  split_RDAforest_spatial <- function(combined_fit, period_indices) {
    retained_global <- which(combined_fit$goodrows)
    prediction_rows <- which(retained_global %in% period_indices)
    result <- combined_fit
    result$predictions.direct <- combined_fit$predictions.direct[
      prediction_rows, , drop = FALSE
    ]
    result$predictions.turnover <- combined_fit$predictions.turnover[
      prediction_rows, , drop = FALSE
    ]
    result$goodrows <- combined_fit$goodrows[period_indices]
    result
  }

  period_n <- nrow(RDAforest_spatial_xy)
  present_indices <- seq_len(period_n)
  mid_indices <- period_n + seq_len(period_n)
  late_indices <- 2 * period_n + seq_len(period_n)
  RDAforest_present_spatial <- split_RDAforest_spatial(
    RDAforest_combined_spatial, present_indices
  )
  RDAforest_mid_spatial <- split_RDAforest_spatial(
    RDAforest_combined_spatial, mid_indices
  )
  RDAforest_late_spatial <- split_RDAforest_spatial(
    RDAforest_combined_spatial, late_indices
  )

  stopifnot(RDAforest_present_spatial$failed.reps == 0)
  stopifnot(RDAforest_mid_spatial$failed.reps == 0)
  stopifnot(RDAforest_late_spatial$failed.reps == 0)

  RDAforest_grid_dx_observed <- median(diff(sort(unique(
    RDAforest_spatial_xy$x
  ))))
  RDAforest_grid_dy_observed <- median(diff(sort(unique(
    RDAforest_spatial_xy$y
  ))))
  RDAforest_grid_anisotropy <-
    (RDAforest_grid_dx_observed * cos(mean(RDAforest_spatial_xy$y) * pi / 180)) /
    RDAforest_grid_dy_observed
  RDAforest_artifact_audit <- data.frame(
    metric = c(
      "grid_rows_terrestrial",
      "unique_longitudes",
      "unique_latitudes",
      "duplicate_coordinates",
      "median_longitude_spacing_degrees",
      "median_latitude_spacing_degrees",
      "east_west_to_north_south_spacing_ratio",
      "present_prediction_retention",
      "mid_prediction_retention",
      "late_prediction_retention"
    ),
    value = c(
      nrow(RDAforest_spatial_xy),
      length(unique(RDAforest_spatial_xy$x)),
      length(unique(RDAforest_spatial_xy$y)),
      sum(duplicated(RDAforest_spatial_xy)),
      RDAforest_grid_dx_observed,
      RDAforest_grid_dy_observed,
      RDAforest_grid_anisotropy,
      mean(RDAforest_present_spatial$goodrows),
      mean(RDAforest_mid_spatial$goodrows),
      mean(RDAforest_late_spatial$goodrows)
    )
  )
  write.csv(
    RDAforest_artifact_audit,
    file.path(RDAforest_spatial_out, "21_RDAforest_spatial_artifact_audit.csv"),
    row.names = FALSE
  )

  set.seed(4040)
  RDAforest_adaptation_q90 <- as.numeric(RDAforest::adapt_scale(
    RDAforest_present_spatial$predictions.direct
  )[2])

  RDAforest_offset_mid <- RDAforest::gen_offset_oj(
    X = RDAforest_present_spatial,
    Y = RDAforest_mid_spatial,
    sx = RDAforest_spatial_xy,
    sy = RDAforest_spatial_xy,
    sc = RDAforest_adaptation_q90
  )
  RDAforest_offset_late <- RDAforest::gen_offset_oj(
    X = RDAforest_present_spatial,
    Y = RDAforest_late_spatial,
    sx = RDAforest_spatial_xy,
    sy = RDAforest_spatial_xy,
    sc = RDAforest_adaptation_q90
  )

  RDAforest_offset_mid$period <- "2041-2070"
  RDAforest_offset_late$period <- "2071-2100"
  RDAforest_offsets <- rbind(RDAforest_offset_mid, RDAforest_offset_late)
  write.csv(
    RDAforest_offsets,
    file.path(RDAforest_spatial_out, "21_RDAforest_original_offset_pixels.csv"),
    row.names = FALSE
  )

  compare_RDAforest_offset <- function(offset_pixels, GF_offset, RDA_offset) {
    offset_points <- terra::vect(
      offset_pixels,
      geom = c("x", "y"),
      crs = terra::crs(GF_offset)
    )
    extracted <- terra::extract(c(GF_offset, RDA_offset), offset_points)
    comparison <- data.frame(
      x = offset_pixels$x,
      y = offset_pixels$y,
      RDAforest = offset_pixels$offset,
      Gradient_Forest = extracted[[2]],
      RDA = extracted[[3]]
    )
    comparison <- comparison[complete.cases(comparison), , drop = FALSE]
    list(
      data = comparison,
      rho_GF = cor(
        comparison$RDAforest, comparison$Gradient_Forest,
        method = "spearman"
      ),
      rho_RDA = cor(
        comparison$RDAforest, comparison$RDA,
        method = "spearman"
      )
    )
  }

  RDAforest_comparison_mid <- compare_RDAforest_offset(
    RDAforest_offset_mid, GF_offset_mid, RDA_offset_mid
  )
  RDAforest_comparison_late <- compare_RDAforest_offset(
    RDAforest_offset_late, GF_offset_late, RDA_offset_late
  )

  RDAforest_offset_summary <- data.frame(
    period = c("2041-2070", "2071-2100"),
    n_terrestrial_cells = c(
      nrow(RDAforest_offset_mid), nrow(RDAforest_offset_late)
    ),
    mean_q90_scaled_offset = c(
      mean(RDAforest_offset_mid$offset),
      mean(RDAforest_offset_late$offset)
    ),
    median_q90_scaled_offset = c(
      median(RDAforest_offset_mid$offset),
      median(RDAforest_offset_late$offset)
    ),
    spearman_with_GF = c(
      RDAforest_comparison_mid$rho_GF,
      RDAforest_comparison_late$rho_GF
    ),
    spearman_with_RDA = c(
      RDAforest_comparison_mid$rho_RDA,
      RDAforest_comparison_late$rho_RDA
    ),
    present_adaptation_q90 = RDAforest_adaptation_q90
  )
  write.csv(
    RDAforest_offset_summary,
    file.path(RDAforest_spatial_out, "21_RDAforest_original_offset_summary.csv"),
    row.names = FALSE
  )

  # Compare within-period percentile ranks because the three offsets use different units.
  # Contrast the RDA rank with the mean rank of GF and RDAforest.
  rank01 <- function(x) {
    if (length(x) <= 1) return(rep(NA_real_, length(x)))
    (rank(x, ties.method = "average") - 1) / (length(x) - 1)
  }

  RDAforest_concordance_mid <- transform(
    RDAforest_comparison_mid$data,
    period = "2041-2070"
  )
  RDAforest_concordance_late <- transform(
    RDAforest_comparison_late$data,
    period = "2071-2100"
  )
  RDAforest_concordance <- rbind(
    RDAforest_concordance_mid,
    RDAforest_concordance_late
  )
  RDAforest_concordance$linear_rank <- ave(
    RDAforest_concordance$RDA,
    RDAforest_concordance$period,
    FUN = rank01
  )
  RDAforest_concordance$GF_rank <- ave(
    RDAforest_concordance$Gradient_Forest,
    RDAforest_concordance$period,
    FUN = rank01
  )
  RDAforest_concordance$RDAforest_rank <- ave(
    RDAforest_concordance$RDAforest,
    RDAforest_concordance$period,
    FUN = rank01
  )
  RDAforest_concordance$nonlinear_rank <- rowMeans(
    RDAforest_concordance[, c("GF_rank", "RDAforest_rank")]
  )
  RDAforest_concordance$consensus_rank <- rowMeans(
    RDAforest_concordance[, c("linear_rank", "nonlinear_rank")]
  )
  RDAforest_concordance$agreement <- 1 - abs(
    RDAforest_concordance$linear_rank -
      RDAforest_concordance$nonlinear_rank
  )
  RDAforest_concordance$interpretation <- "concordant intermediate"
  RDAforest_concordance$interpretation[
    RDAforest_concordance$linear_rank >= .67 &
      RDAforest_concordance$nonlinear_rank >= .67
  ] <- "concordant high"
  RDAforest_concordance$interpretation[
    RDAforest_concordance$linear_rank <= .33 &
      RDAforest_concordance$nonlinear_rank <= .33
  ] <- "concordant low"
  RDAforest_concordance$interpretation[
    RDAforest_concordance$linear_rank -
      RDAforest_concordance$nonlinear_rank > .25
  ] <- "linear higher"
  RDAforest_concordance$interpretation[
    RDAforest_concordance$nonlinear_rank -
      RDAforest_concordance$linear_rank > .25
  ] <- "nonlinear higher"
  RDAforest_concordance$interpretation <- factor(
    RDAforest_concordance$interpretation,
    levels = c(
      "concordant low", "concordant intermediate", "concordant high",
      "linear higher", "nonlinear higher"
    )
  )
  write.csv(
    RDAforest_concordance,
    file.path(
      RDAforest_spatial_out,
      "21_linear_nonlinear_offset_concordance_pixels.csv"
    ),
    row.names = FALSE
  )

  concordance_late_raster <- terra::rast(
    RDAforest_concordance[
      RDAforest_concordance$period == "2071-2100",
      c(
        "x", "y", "linear_rank", "GF_rank", "RDAforest_rank",
        "nonlinear_rank", "consensus_rank", "agreement"
      )
    ],
    type = "xyz",
    crs = terra::crs(GF_offset_late)
  )
  writeRaster(
    concordance_late_raster,
    file.path(
      RDAforest_spatial_out,
      "21_linear_nonlinear_offset_concordance_late.tif"
    ),
    overwrite = TRUE
  )

  present_good <- RDAforest_present_spatial$goodrows
  present_xy <- RDAforest_spatial_xy[present_good, , drop = FALSE]
  present_environment <- RDAforest_present_X[present_good, , drop = FALSE]

  # Obtain the predicted genomic-space ordination for the display panels.
  grDevices::pdf(NULL)
  RDAforest_adaptation_unclustered <- RDAforest::plot_adaptation(
    RDAforest_present_spatial$predictions.direct,
    present_environment,
    present_xy,
    main = "RDAforest present-day adaptation",
    rangeExp = 1.5, scal = 10, jitscale = .05,
    color.scheme = "001", lighten = .65
  )
  dev.off()
  RDAforest_adaptation_panels <- adaptation_ggplot_panels(
    RDAforest_adaptation_unclustered,
    present_environment
  )
  save_ggplot_row(
    RDAforest_adaptation_panels,
    file.path(figures, "22_RDAforest_adaptation_unclustered.png"),
    width = 2400, height = 1200, res = 220
  )

  grDevices::pdf(NULL)
  RDAforest_neighborhoods_direct <- RDAforest::plot_adaptation(
    RDAforest_present_spatial$predictions.direct,
    present_environment,
    present_xy,
    main = "Direct-prediction neighborhoods",
    rangeExp = 1.5, scal = 10, jitscale = .05,
    color.scheme = "001", lighten = .65,
    nclust = 12, cluster.merge.level = .25
  )
  dev.off()
  RDAforest_direct_panels <- adaptation_ggplot_panels(
    RDAforest_neighborhoods_direct,
    present_environment
  )
  save_ggplot_row(
    RDAforest_direct_panels,
    file.path(figures, "23_RDAforest_neighborhoods_direct.png"),
    width = 3300, height = 1100, res = 200
  )

  grDevices::pdf(NULL)
  RDAforest_neighborhoods_turnover <- RDAforest::plot_adaptation(
    RDAforest_present_spatial$predictions.direct,
    present_environment,
    present_xy,
    cluster.guide = RDAforest_present_spatial$predictions.turnover,
    main = "Turnover-guided neighborhoods",
    rangeExp = 1.5, scal = 10, jitscale = .05,
    color.scheme = "001", lighten = .65,
    nclust = 12, cluster.merge.level = .25
  )
  dev.off()
  RDAforest_turnover_panels <- adaptation_ggplot_panels(
    RDAforest_neighborhoods_turnover,
    present_environment
  )
  save_ggplot_row(
    RDAforest_turnover_panels,
    file.path(figures, "24_RDAforest_neighborhoods_turnover.png"),
    width = 3300, height = 1100, res = 200
  )

  # Fit the display PCA to present-day direct predictions and project both future periods.
  RDAforest_prediction_pca <- prcomp(
    RDAforest_present_spatial$predictions.direct,
    center = TRUE, scale. = TRUE
  )
  project_RDAforest_prediction <- function(spatial_fit) {
    predict(
      RDAforest_prediction_pca,
      newdata = spatial_fit$predictions.direct
    )[, 1:3, drop = FALSE]
  }
  RDAforest_present_scores <- project_RDAforest_prediction(
    RDAforest_present_spatial
  )
  RDAforest_mid_scores <- project_RDAforest_prediction(RDAforest_mid_spatial)
  RDAforest_late_scores <- project_RDAforest_prediction(RDAforest_late_spatial)
  # Use present-day empirical percentiles for the three axes in every period.
  RDAforest_axis_ECDF <- lapply(
    seq_len(ncol(RDAforest_present_scores)),
    function(j) stats::ecdf(RDAforest_present_scores[, j])
  )
  RDAforest_IDV_colour <- function(scores) {
    scores01 <- vapply(
      seq_len(ncol(scores)),
      function(j) RDAforest_axis_ECDF[[j]](scores[, j]),
      numeric(nrow(scores))
    )
    idv_false_colour(
      scores01[, 1], scores01[, 2], scores01[, 3]
    )
  }

  RDAforest_colour_map_data <- rbind(
    transform(
      RDAforest_spatial_xy[
        RDAforest_present_spatial$goodrows, , drop = FALSE
      ],
      brand_colour = RDAforest_IDV_colour(RDAforest_present_scores),
      period = "Present"
    ),
    transform(
      RDAforest_spatial_xy[
        RDAforest_mid_spatial$goodrows, , drop = FALSE
      ],
      brand_colour = RDAforest_IDV_colour(RDAforest_mid_scores),
      period = "2041-2070"
    ),
    transform(
      RDAforest_spatial_xy[
        RDAforest_late_spatial$goodrows, , drop = FALSE
      ],
      brand_colour = RDAforest_IDV_colour(RDAforest_late_scores),
      period = "2071-2100"
    )
  )
  RDAforest_colour_map_data$period <- factor(
    RDAforest_colour_map_data$period,
    levels = c("Present", "2041-2070", "2071-2100")
  )
  RDAforest_colour_map_data$axis_key <- NA_character_
  RDAforest_axis_colours <- c(
    "Axis 1 dominant" = "#6A3D9A",
    "Axis 2 dominant" = "#E69F00",
    "Axis 3 dominant" = "#00796B"
  )
  RDAforest_grid_dx <- median(diff(sort(unique(RDAforest_colour_map_data$x))))
  RDAforest_grid_dy <- median(diff(sort(unique(RDAforest_colour_map_data$y))))
  plot_RDAforest_colour_map <- ggplot(
    RDAforest_colour_map_data,
    aes(x, y, fill = brand_colour)
  ) +
    geom_polygon(data = ameiva_world_map, aes(long, lat, group = group), inherit.aes = FALSE, fill = idv_cream, colour = "grey78", linewidth = .2) +
    geom_tile(
      width = RDAforest_grid_dx * 1.015,
      height = RDAforest_grid_dy * 1.015,
      linewidth = 0
    ) +
    scale_fill_identity(
      name = "Dominant genomic axis",
      breaks = unname(RDAforest_axis_colours),
      labels = names(RDAforest_axis_colours),
      guide = guide_legend(override.aes = list(colour = NA))
    ) +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    theme_ameiva() +
    theme(panel.grid = element_blank()) +
    labs(
      title = "RDAforest adaptive landscape through time",
      subtitle = paste0(
        "ECDF-transformed axes: purple = axis 1, orange = axis 2, ",
        "blue-green = axis 3; mixtures indicate joint responses"
      ),
      x = "Longitude", y = "Latitude"
    )
  ggsave(
    file.path(figures, "25_RDAforest_present_future_RGB.png"),
    plot_RDAforest_colour_map, width = 13, height = 5,
    dpi = 400, bg = "white"
  )
  ggsave(
    file.path(figures, "25_RDAforest_present_future_RGB.tiff"),
    plot_RDAforest_colour_map, width = 13, height = 5,
    dpi = 600, bg = "white", compression = "lzw"
  )
  ggsave(
    file.path(figures, "25_RDAforest_present_future_RGB.pdf"),
    plot_RDAforest_colour_map, width = 13, height = 5,
    bg = "white", device = cairo_pdf
  )

  # Use square points to represent the grid despite floating-point coordinate differences.
  gc()
  plot_RDAforest_offset <- ggplot(
    RDAforest_offsets,
    aes(x, y, colour = offset)
  ) +
    geom_point(shape = 15, size = .62) +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap() +
    scale_colour_gradientn(colours = idv_sequential) +
    theme_ameiva() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "RDAforest genomic offset",
      subtitle = "gen_offset_oj distance scaled by present-day q90",
      colour = "q90-scaled\noffset"
    )
  ggsave(file.path(figures, "26_RDAforest_original_genomic_offset.png"),
         plot_RDAforest_offset, width = 11, height = 5.5,
         dpi = 400, bg = "white")

  RDAforest_offset_comparison_plot <- rbind(
    transform(RDAforest_comparison_mid$data, period = "2041-2070"),
    transform(RDAforest_comparison_late$data, period = "2071-2100")
  )
  plot_RDAforest_comparison <- ggplot(
    RDAforest_offset_comparison_plot,
    aes(RDAforest, Gradient_Forest)
  ) +
    geom_point(alpha = .10, size = .45, colour = idv_purple) +
    stat_summary_bin(
      fun = median, geom = "line", bins = 35,
      colour = idv_orange, linewidth = 1
    ) +
    facet_wrap(~period, scales = "free") +
    theme_ameiva() +
    labs(
      title = "Agreement between independently constructed genomic offsets",
      subtitle = sprintf(
        "Spearman rho: %.2f (2041-2070) and %.2f (2071-2100)",
        RDAforest_comparison_mid$rho_GF,
        RDAforest_comparison_late$rho_GF
      ),
      x = "RDAforest offset (present q90 units)",
      y = "Gradient Forest transformed-space offset"
    )
  ggsave(file.path(figures, "27_RDAforest_GF_offset_comparison.png"),
         plot_RDAforest_comparison, width = 10, height = 5.5,
         dpi = 400, bg = "white")

  # Map relative offsets on the shared valid domain, not raw distances across methods.
  concordance_late <- subset(
    RDAforest_concordance, period == "2071-2100"
  )
  concordance_colours <- c(
    "concordant low" = idv_purple_dark,
    "concordant intermediate" = idv_teal,
    "concordant high" = idv_orange_dark,
    "linear higher" = idv_purple,
    "nonlinear higher" = idv_orange
  )
  concordance_tile_width <- 1.02 * median(diff(sort(unique(concordance_late$x))))
  concordance_tile_height <- 1.02 * median(diff(sort(unique(concordance_late$y))))

  plot_linear_rank <- ggplot(
    concordance_late, aes(x, y, fill = linear_rank)
  ) +
    geom_tile(width = concordance_tile_width, height = concordance_tile_height) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    scale_fill_gradientn(colours = idv_sequential, limits = c(0, 1)) +
    theme_ameiva() + theme(panel.grid = element_blank()) +
    labs(
      title = "Linear climate association",
      subtitle = "Partial-RDA offset rank, 2071-2100",
      x = "Longitude", y = "Latitude", fill = "Continental\nrank"
    )

  plot_GF_rank <- ggplot(
    concordance_late, aes(x, y, fill = GF_rank)
  ) +
    geom_tile(width = concordance_tile_width, height = concordance_tile_height) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    scale_fill_gradientn(colours = idv_sequential, limits = c(0, 1)) +
    theme_ameiva() + theme(panel.grid = element_blank()) +
    labs(
      title = "Candidate nonlinear association",
      subtitle = "Gradient Forest offset rank, 2071-2100",
      x = "Longitude", y = "Latitude", fill = "Continental\nrank"
    )

  plot_RDAforest_rank <- ggplot(
    concordance_late, aes(x, y, fill = RDAforest_rank)
  ) +
    geom_tile(width = concordance_tile_width, height = concordance_tile_height) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    scale_fill_gradientn(colours = idv_sequential, limits = c(0, 1)) +
    theme_ameiva() + theme(panel.grid = element_blank()) +
    labs(
      title = "Genome-wide nonlinear association",
      subtitle = "RDAforest offset rank, 2071-2100",
      x = "Longitude", y = "Latitude", fill = "Continental\nrank"
    )

  plot_linear_nonlinear_concordance <- ggplot(
    concordance_late, aes(x, y, fill = interpretation)
  ) +
    geom_tile(width = concordance_tile_width, height = concordance_tile_height) +
    coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15)) +
    scale_fill_manual(values = concordance_colours, drop = FALSE) +
    theme_ameiva() + theme(panel.grid = element_blank()) +
    labs(
      title = "Linear–nonlinear concordance",
      subtitle = "Agreement of relative displacement, 2071-2100",
      x = "Longitude", y = "Latitude", fill = "Spatial reading"
    )

  plot_offset_concordance <- patchwork::wrap_plots(
    plot_linear_rank,
    plot_GF_rank,
    plot_RDAforest_rank,
    plot_linear_nonlinear_concordance,
    ncol = 2
  ) + patchwork::plot_annotation(tag_levels = "a")
  ggsave(
    file.path(figures, "35_linear_nonlinear_offset_concordance.png"),
    plot_offset_concordance, width = 14, height = 10,
    dpi = 400, bg = "white"
  )

  plot_RDAforest_domain <- ggplot(
    subset(RDAforest_domain_audit, period != "Present"),
    aes(x, y, colour = status)
  ) +
    geom_point(shape = 15, size = .55) +
    facet_wrap(~period, nrow = 1) +
    coord_quickmap() +
    scale_colour_manual(values = c(
      inside = idv_purple,
      clamped = idv_orange,
      excluded = idv_grey
    )) +
    theme_ameiva() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "RDAforest projection-domain audit",
      subtitle = "Purple: inside training range; orange: admitted and clamped; grey: excluded",
      colour = "Projection status"
    )
  ggsave(file.path(figures, "28_RDAforest_projection_domain_audit.png"),
         plot_RDAforest_domain, width = 11, height = 5.5,
         dpi = 400, bg = "white")
}


# 22. Assemble the four composite figures from the core workflow ----

# Export these composite panels and retain the individual diagnostic figures.
if (
  RUN_GRADIENT_FOREST && RUN_RDAFOREST && RUN_SPATIAL_PROJECTION &&
  RUN_RDAFOREST_SPATIAL
) {
  main_figure_1 <- patchwork::wrap_plots(
    plot_population_map,
    plot_pca,
    plot_snmf_ancestry,
    ncol = 2,
    heights = c(1, .72),
    design = "AB\nCC"
  ) +
    patchwork::plot_annotation(tag_levels = "a")

  main_figure_2 <- patchwork::wrap_plots(
    plot_partial_RDA,
    plot_mosaic_PCA,
    plot_mosaic_map,
    design = "AABB\nAABB\nCCCC"
  ) +
    patchwork::plot_annotation(tag_levels = "a")

  main_figure_3 <- patchwork::wrap_plots(
    plot_GF_importance,
    plot_importance_uncertainty,
    plot_candidate_function,
    design = "AABB\nAABB\nCCCC"
  ) +
    patchwork::plot_annotation(tag_levels = "a")

  main_figure_4 <- plot_offset_concordance

  main_figures <- list(
    "Main_Figure_1_structure" = main_figure_1,
    "Main_Figure_2_linear_adaptive_mosaic" = main_figure_2,
    "Main_Figure_3_polygenic_functional_architecture" = main_figure_3,
    "Main_Figure_4_landscape_and_offset" = main_figure_4
  )
  main_figure_dimensions <- data.frame(
    figure = names(main_figures),
    width = c(13, 16, 16, 14),
    height = c(10, 11, 10, 10)
  )

  for (i in seq_along(main_figures)) {
    figure_name <- names(main_figures)[i]
    figure_width <- main_figure_dimensions$width[i]
    figure_height <- main_figure_dimensions$height[i]
    ggsave(
      file.path(figures, paste0(figure_name, ".png")),
      main_figures[[i]], width = figure_width, height = figure_height,
      dpi = 400, bg = "white", limitsize = FALSE
    )
    ggsave(
      file.path(figures, paste0(figure_name, ".tiff")),
      main_figures[[i]], width = figure_width, height = figure_height,
      dpi = 400, bg = "white", compression = "lzw", limitsize = FALSE
    )
  }

  write.csv(
    data.frame(
      figure = names(main_figures),
      panels = c(3, 3, 3, 4),
      role = c(
        "continental sampling, genomic PCA and continuous K=4 ancestry",
        "partial RDA with allelic response and the consensus candidate mosaic",
        "candidate-focused GF, genome-wide RDAforest and functional context",
        "2071-2100 RDA, GF and RDAforest ranks plus linear-nonlinear concordance"
      )
    ),
    file.path(out, "00_MAIN_FIGURE_REGISTRY.csv"),
    row.names = FALSE
  )
}


# 23. Export analysis summaries, figure index and session information ----

summary_final <- data.frame(
  metric = c(
    "individuals",
    "SNPs_main_panel",
    "RAD_loci_represented",
    "sNMF_bestK",
    "environmental_predictors",
    "RDA_Josue_candidates",
    "RDA_PCA_geo_candidates",
    "RDA_sNMF_geo_candidates",
    "RDA_landscape_consensus",
    "RDA_landscape_strict_3SD_consensus",
    "LFMM_candidates",
    "GEA_consensus_candidates",
    "Gradient_Forest_loci",
    "adaptive_mosaic_K",
    "RDA_landscape_K",
    "GF_landscape_K",
    "RDA_GF_spearman_mid_century",
    "RDA_GF_spearman_late_century",
    "RDAforest_selected_predictors",
    "RDAforest_IBD_protest_r",
    "RDAforest_IBD_permutation_P",
    "RDAforest_residual_IBD_mantel_spearman",
    "RDAforest_residual_IBD_permutation_P",
    "RDAforest_median_total_importance",
    "RDAforest_offset_mean_mid_century",
    "RDAforest_offset_mean_late_century",
    "RDAforest_GF_spearman_mid_century",
    "RDAforest_GF_spearman_late_century",
    "RDAforest_RDA_spearman_mid_century",
    "RDAforest_RDA_spearman_late_century"
  ),
  value = c(
    nrow(geno_complete),
    ncol(geno_complete),
    length(unique(rad_locus)),
    bestK,
    length(variables_final),
    length(candidates_RDA),
    length(candidates_RDA_pca_geo),
    length(candidates_RDA_snmf_geo),
    length(candidates_RDA_material),
    length(candidates_RDA_material_strict),
    length(candidates_LFMM),
    length(candidates_consensus),
    length(GF_loci),
    best_mosaic_k,
    if (RUN_SPATIAL_PROJECTION) RDA_landscape$selected_k else NA,
    if (RUN_SPATIAL_PROJECTION) GF_landscape$selected_k else NA,
    if (RUN_SPATIAL_PROJECTION) offset_summary$RDA_GF_spearman[1] else NA,
    if (RUN_SPATIAL_PROJECTION) offset_summary$RDA_GF_spearman[2] else NA,
    if (RUN_RDAFOREST) length(selected_RDAforest) else NA,
    if (RUN_RDAFOREST) RDAforest_IBD_test$t0 else NA,
    if (RUN_RDAFOREST) RDAforest_IBD_test$signif else NA,
    if (RUN_RDAFOREST) unname(RDAforest_residual_IBD_test$statistic) else NA,
    if (RUN_RDAFOREST) RDAforest_residual_IBD_test$signif else NA,
    if (RUN_RDAFOREST) median(rdaforest_fit$sum.importances) else NA,
    if (RUN_RDAFOREST_SPATIAL) RDAforest_offset_summary$mean_q90_scaled_offset[1] else NA,
    if (RUN_RDAFOREST_SPATIAL) RDAforest_offset_summary$mean_q90_scaled_offset[2] else NA,
    if (RUN_RDAFOREST_SPATIAL) RDAforest_offset_summary$spearman_with_GF[1] else NA,
    if (RUN_RDAFOREST_SPATIAL) RDAforest_offset_summary$spearman_with_GF[2] else NA,
    if (RUN_RDAFOREST_SPATIAL) RDAforest_offset_summary$spearman_with_RDA[1] else NA,
    if (RUN_RDAFOREST_SPATIAL) RDAforest_offset_summary$spearman_with_RDA[2] else NA
  )
)

write.csv(summary_final, file.path(out, "00_FINAL_CHECKPOINT.csv"), row.names = FALSE)

plot_registry <- data.frame(
  plot = c(
    "01_sNMF_cross_entropy.png",
    "02_PCA_sNMF_structure.png",
    "03_Gradient_Forest_importance.png",
    "04_RDAforest_importance.png",
    "05_adaptive_mosaic_genotype_PCA.png",
    "06_adaptive_mosaic_map.png",
    "07_GF_present_RGB.png",
    "08_GF_future_RGB.png",
    "09_GF_genomic_offset.png",
    "10_RDA_present_landscape.png",
    "11_RDA_genomic_offset.png",
    "12_RDA_GF_offset_disagreement.png",
    "13_landscape_kmeans_RDA_GF.png",
    "14_RDAforest_IBD_diagnostic.png",
    "15_RDAforest_axes_predictability.png",
    "16_RDAforest_mtry_delta.png",
    "17_RDAforest_mtry_positive.png",
    "18_RDAforest_importance_uncertainty.png",
    "19_RDAforest_exploratory_turnover.png",
    "20_RDAforest_jackknife_turnover.png",
    "21_RDAforest_sample_adaptation.png",
    "22_RDAforest_adaptation_unclustered.png",
    "23_RDAforest_neighborhoods_direct.png",
    "24_RDAforest_neighborhoods_turnover.png",
    "25_RDAforest_present_future_RGB.png",
    "26_RDAforest_original_genomic_offset.png",
    "27_RDAforest_GF_offset_comparison.png",
    "28_RDAforest_projection_domain_audit.png",
    "29_GF_RDA_offset_domain_audit.png",
    "30_climate_change_input_artifact_audit.png",
    "31_partial_RDA_genotype_response.png",
    "32_population_sampling_map.png",
    "33_sNMF_ancestry_barplot.png",
    "34_candidate_functional_context.png"
  ),
  script_section = c(
    6, 7, 15, 16, 17, 17, 18, 18, 19, 18, 19, 19, 20,
    rep(16, 8), rep(21, 7), 20, 20, 11, 7, 7, 14
  ),
  biological_role = c(
    "select demographic complexity",
    "describe population structure",
    "identify nonlinear environmental turnover",
    "identify predictors of total nonlinear genomic covariance",
    "cluster adaptive multilocus genotypes",
    "map the adaptive mosaic within Ameiva ameiva",
    "map present nonlinear adaptive turnover in project-constrained false colour",
    "project nonlinear adaptive turnover into future periods",
    "map transformed-space Gradient Forest genomic offset",
    "map the present linear adaptive landscape along RDA1",
    "map eigenvalue-weighted RDA genomic offset",
    "locate disagreement between linear and nonlinear offsets",
    "compare independently inferred RDA and Gradient Forest landscape regions",
    "diagnose isolation by distance before RDAforest conditioning",
    "show residual genomic inertia and out-of-bag predictability of genomic axes",
    "show predictor importance change between low and high mtry",
    "document which predictors pass the RDAforest mtry criterion",
    "quantify ordination-jackknife uncertainty in environmental importance",
    "show exploratory cumulative turnover curves from makeGF",
    "show eigenvalue-scaled turnover curves from ordinationJackknife",
    "map sample-level direct predictions in genomic-PC colour space",
    "map continuous present-day RDAforest adaptation",
    "map adaptive neighborhoods clustered from direct predictions",
    "map adaptive neighborhoods guided by turnover curves",
    "compare present and future RDAforest prediction colours",
    "map q90-scaled genomic offset from the original gen_offset_oj function",
    "compare the independently constructed RDAforest and Gradient Forest offsets",
    "audit clamping and exclusion outside the RDAforest training domain",
    "audit the four-variable environmental domain retained for GF and RDA offsets",
    "show which map patterns are already present in future-minus-present climate inputs",
    "show the partial RDA with the allelic matrix explicitly identified as response",
    "map sampling and dominant continuous ancestry components",
    "show individual sNMF K=4 ancestry coefficients",
    "summarize the genomic functional context of RDA-and-LFMM candidates"
  )
)

write.csv(plot_registry, file.path(out, "00_PLOT_REGISTRY.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "00_sessionInfo.txt"))

print(summary_final)
