
# 1. Read signed LFMM2 effects, genotypes, environment and spatial folds ----

library(ggplot2)

base <- Sys.getenv("AMEIVA_DEVELOPMENT_DIR")
results <- Sys.getenv("AMEIVA_RESULTS_DIR")
project <- Sys.getenv("AMEIVA_PROJECT_DIR")
environment_file <- Sys.getenv("AMEIVA_ENV_FILE")
stopifnot(nzchar(base), nzchar(results), nzchar(project),
          nzchar(environment_file), dir.exists(base), dir.exists(results),
          file.exists(environment_file))

# Confirm the reconstructed LFMM2 effects reproduce the saved test statistics.
check <- read.csv(file.path(base, "lfmm2_recovery/01_refit_vs_original_test_check.csv"))
stopifnot(nrow(check) == 4L, all(check$n_SNPs == 14134L),
          all(check$max_abs_z_difference < 1e-8),
          all(check$max_abs_p_difference < 1e-8),
          all(check$coefficient_z_sign_agreement == 1))
groups <- read.csv(file.path(base, "14_LFMM2_signed_profile_groups_v2.csv"))
effects <- read.csv(file.path(base, "15_LFMM2_consensus_effect_profiles_v2.csv"))
site_scores <- read.csv(file.path(base, "18_LFMM2_group_locality_scores_v2.csv"),
                        check.names = FALSE)
candidate_scores <- read.csv(file.path(base, "05_individual_multilocus_scores.csv"),
                             check.names = FALSE)
fold_refit <- read.csv(file.path(base, "22_LFMM2_spatial_group_stability_v3.csv"))
folds <- read.csv(file.path(project,
  "results/landscape/checkerboard4_folds/02_checkerboard4_by_individual.csv"))
G <- readRDS(file.path(results, "03_genotype_matrix_14134.rds"))
geno <- as.matrix(read.table(file.path(results,
  "03_Ameiva_14134.lfmm_imputed.lfmm"), header = FALSE, check.names = FALSE))
rownames(geno) <- rownames(G$geno)
colnames(geno) <- colnames(G$geno)
PCA <- read.csv(file.path(results, "04_PCA_14134_scores.csv"))
Variables <- read.csv(environment_file, check.names = FALSE)
# Match all individual-level tables to the genotype matrix by identifier.
PCA <- PCA[match(rownames(geno), PCA$InternalID), ]
Variables <- Variables[match(rownames(geno), Variables$InternalID), ]
folds <- folds[match(rownames(geno), folds$InternalID), ]
candidate_scores <- candidate_scores[match(rownames(geno),
                                            candidate_scores$InternalID), ]
RDA_ids <- unique(read.csv(file.path(results,
  "08_RDA_candidates_landscape_consensus.csv"))$stacks_id)
LFMM_ids <- unique(read.csv(file.path(results,
  "09_LFMM2_candidates_FDR05.csv"))$stacks_id)
stopifnot(identical(dim(geno), c(245L, 14134L)),
          identical(rownames(geno), PCA$InternalID),
          identical(rownames(geno), Variables$InternalID),
          identical(rownames(geno), folds$InternalID),
          identical(rownames(geno), candidate_scores$InternalID),
          all(is.finite(geno)), all(geno %in% 0:2),
          length(intersect(RDA_ids, LFMM_ids)) == 164L,
          setequal(intersect(RDA_ids, LFMM_ids), groups$stacks_id),
          !anyDuplicated(groups$stacks_id),
          all(groups$stacks_id %in% colnames(geno)),
          length(unique(candidate_scores$locality)) == 185L,
          setequal(unique(groups$group), c("BIO2 +", "BIO2 -", "BIO3 +",
                                        "BIO3 -", "BIO13 +", "BIO13 -",
                                        "BIO14 +", "BIO14 -")))

# 2. Summarise SNP and RAD-locus counts and signed LFMM2 effects ----

count <- as.data.frame(table(groups$group), stringsAsFactors = FALSE)
names(count) <- c("group", "n_SNPs")
count$n_RAD_loci <- vapply(count$group, function(g)
  length(unique(groups$rad_locus[groups$group == g])), integer(1))
count <- count[order(count$n_SNPs, decreasing = TRUE), ]
write.csv(count, file.path(base, "30_v4_group_counts.csv"), row.names = FALSE)

heat <- aggregate(standardized_beta ~ group + predictor, effects, mean)
heat$group_label <- paste0(heat$group, "  ",
  count$n_SNPs[match(heat$group, count$group)], "/",
  count$n_RAD_loci[match(heat$group, count$group)])
heat$group_label <- factor(heat$group_label,
  levels = paste0(count$group, "  ", count$n_SNPs, "/", count$n_RAD_loci))
heat$predictor <- factor(toupper(heat$predictor),
                         levels = c("BIO2", "BIO3", "BIO13", "BIO14"))
write.csv(heat, file.path(base, "31_v4_group_mean_signed_effects.csv"),
          row.names = FALSE)

# 3. Estimate group mean ALT dosage by individual and locality ----

climate_names <- c(BIO2 = "Diurnal temp. range (°C)",
                   BIO3 = "Isothermality (ratio)",
                   BIO13 = "Wettest-month precip. (mm)",
                   BIO14 = "Driest-month precip. (mm)")
frequency <- data.frame()
individual <- data.frame()
# Dividing mean diploid dosage by two expresses each group on an allele-frequency scale.
for (g in count$group) {
  loci <- groups$stacks_id[groups$group == g]
  y <- rowMeans(geno[, loci, drop = FALSE]) / 2
  predictor <- sub(" [+-]$", "", g)
  x <- Variables[[paste0("X.CHELSA_", tolower(predictor),
                          "_1981.2010_SA")]]
  stopifnot(all(is.finite(y)), all(y >= 0 & y <= 1), all(is.finite(x)))
  individual <- rbind(individual,
    data.frame(group = g, InternalID = rownames(geno),
               climate_original = x, mean_ALT_frequency = y))
  for (loc in unique(candidate_scores$locality)) {
    members <- which(candidate_scores$locality == loc)
    stopifnot(length(unique(Variables$long[members])) == 1L,
              length(unique(Variables$lat[members])) == 1L)
    frequency <- rbind(frequency,
      data.frame(group = g, predictor = predictor, locality = loc,
                 longitude = Variables$long[members[1]],
                 latitude = Variables$lat[members[1]],
                 POP = Variables$POP[members[1]],
                 n_individuals = length(members),
                 climate_original = mean(x[members]),
                 mean_ALT_frequency = mean(y[members])))
  }
}
stopifnot(nrow(frequency) == 8L * 185L)
write.csv(frequency, file.path(base, "32_v4_all_group_locality_ALT_frequency.csv"),
          row.names = FALSE)

# 4. Summarise locality frequencies across climatic deciles ----

bins <- data.frame()
for (g in count$group) {
  f <- frequency[frequency$group == g, ]
  breaks <- unique(as.numeric(quantile(f$climate_original,
                                       probs = seq(0, 1, .1))))
  f$climate_bin <- cut(f$climate_original, breaks = breaks,
                       include.lowest = TRUE, labels = FALSE)
  stopifnot(all(!is.na(f$climate_bin)))
  for (b in sort(unique(f$climate_bin))) {
    part <- f[f$climate_bin == b, ]
    bins <- rbind(bins,
      data.frame(group = g, predictor = part$predictor[1], bin = b,
                 mean_climate = weighted.mean(part$climate_original,
                                               part$n_individuals),
                 mean_ALT_frequency = weighted.mean(part$mean_ALT_frequency,
                                                    part$n_individuals),
                 n_localities = nrow(part),
                 n_individuals = sum(part$n_individuals)))
  }
}
write.csv(bins, file.path(base, "33_v4_observed_frequency_by_climate_decile.csv"),
          row.names = FALSE)

# 5. Compare conditional associations, linkage sensitivity and spatial holdouts ----

# Fractional-response quasi-binomial GLMs model group frequencies on the 0–1 scale.
# Each individual has equal weight; SNPs are not independent binomial trials.
conditional <- data.frame()
held_out <- data.frame()
rad_check <- data.frame()
for (g in count$group) {
  d <- individual[individual$group == g, ]
  d <- d[match(rownames(geno), d$InternalID), ]
  predictor <- sub(" [+-]$", "", g)
  x <- d$climate_original
  loci <- groups$stacks_id[groups$group == g]
  first_per_RAD <- loci[!duplicated(groups$rad_locus[match(loci,
                                                            groups$stacks_id)])]
  y_one <- rowMeans(geno[, first_per_RAD, drop = FALSE]) / 2
  covariates <- data.frame(y = d$mean_ALT_frequency,
    climate = as.numeric(scale(x)),
    PC1 = PCA$PC1, PC2 = PCA$PC2, PC3 = PCA$PC3,
    lon = Variables$long, lat = Variables$lat)
  fit <- glm(y ~ climate + PC1 + PC2 + PC3 + lon + lat,
             data = covariates, family = quasibinomial(link = "logit"))
  pruned <- covariates
  pruned$y <- y_one
  fit_one <- glm(y ~ climate + PC1 + PC2 + PC3 + lon + lat,
                 data = pruned, family = quasibinomial(link = "logit"))
  stopifnot(isTRUE(fit$converged), isTRUE(fit_one$converged))
  marginal <- lm(mean_ALT_frequency ~ climate_original,
                 data = frequency[frequency$group == g, ],
                 weights = n_individuals)
  conditional <- rbind(conditional,
    data.frame(group = g, predictor = predictor,
      n_SNPs = length(loci), n_RAD_loci = length(first_per_RAD),
      marginal_locality_slope_per_original_unit =
        unname(coef(marginal)["climate_original"]),
      conditional_logit_slope_per_climate_SD =
        unname(coef(fit)["climate"]),
      conditional_direction = sign(coef(fit)["climate"]),
      conditional_slope_agrees_with_LFMM2 =
        sign(coef(fit)["climate"]) == ifelse(grepl("\\+$", g), 1, -1)))
  rad_check <- rbind(rad_check,
    data.frame(group = g, n_SNPs = length(loci),
      n_RAD_loci = length(first_per_RAD),
      correlation_frequency_one_per_RAD = cor(d$mean_ALT_frequency, y_one),
      all_SNP_conditional_logit_slope = unname(coef(fit)["climate"]),
      one_per_RAD_conditional_logit_slope =
        unname(coef(fit_one)["climate"]),
      same_conditional_direction =
        sign(coef(fit)["climate"]) == sign(coef(fit_one)["climate"])))
  # Withhold each checkerboard fold and standardise climate using training individuals.
  # Candidate groups were defined before this diagnostic, so this tests prediction
  # of group frequencies rather than independent candidate discovery.
  for (held in 1:4) {
    train <- folds$fold != held
    data <- covariates
    data$climate <- (x - mean(x[train])) / sd(x[train])
    no_climate <- glm(y ~ PC1 + PC2 + PC3 + lon + lat,
                      data = data[train, ],
                      family = quasibinomial(link = "logit"))
    climate_fit <- glm(y ~ climate + PC1 + PC2 + PC3 + lon + lat,
                       data = data[train, ],
                       family = quasibinomial(link = "logit"))
    stopifnot(isTRUE(no_climate$converged), isTRUE(climate_fit$converged))
    p0 <- predict(no_climate, newdata = data[!train, ], type = "response")
    p1 <- predict(climate_fit, newdata = data[!train, ], type = "response")
    stopifnot(all(is.finite(p0)), all(is.finite(p1)),
              all(p0 >= 0 & p0 <= 1), all(p1 >= 0 & p1 <= 1))
    rmse0 <- sqrt(mean((data$y[!train] - p0)^2))
    rmse1 <- sqrt(mean((data$y[!train] - p1)^2))
    held_out <- rbind(held_out,
      data.frame(group = g, held_out_fold = held,
        n_train = sum(train), n_test = sum(!train),
        training_logit_slope = unname(coef(climate_fit)["climate"]),
        training_direction_matches_LFMM2 =
          sign(coef(climate_fit)["climate"]) ==
            ifelse(grepl("\\+$", g), 1, -1),
        RMSE_nuisance_only = rmse0, RMSE_with_climate = rmse1,
        delta_RMSE = rmse0 - rmse1))
  }
}
write.csv(conditional, file.path(base, "34_v4_marginal_and_conditional_diagnostics.csv"),
          row.names = FALSE)
write.csv(held_out, file.path(base, "35_v4_spatial_holdout_climate_contribution.csv"),
          row.names = FALSE)
write.csv(rad_check, file.path(base, "36_v4_one_SNP_per_RAD_sensitivity.csv"),
          row.names = FALSE)

stability <- aggregate(cbind(group_agrees, top_sign_agrees) ~ original_group,
                       fold_refit, mean)
names(stability)[1] <- "group"
write.csv(stability, file.path(base, "37_v4_LFMM2_group_fold_stability.csv"),
          row.names = FALSE)

# 6. Standardise locality scores and select groups for maps ----

map_data <- do.call(rbind, lapply(count$group, function(g)
  data.frame(site_scores[, c("longitude", "latitude", "n_individuals")],
             group = g, raw_locality_score = site_scores[[g]])))
map_data$display_z <- ave(map_data$raw_locality_score, map_data$group,
                          FUN = function(z) as.numeric(scale(z)))
display_limit <- max(abs(quantile(map_data$display_z, c(.02, .98))))
write.csv(data.frame(scale = "within-group locality z-score", centre = 0,
                     symmetric_display_limit = display_limit,
                     n_saturated = sum(abs(map_data$display_z) > display_limit),
                     n_locality_group_points = nrow(map_data)),
          file.path(base, "38_v4_map_display_scale.csv"), row.names = FALSE)

# Select the two largest groups with at least 25 SNPs and positive predictive gain in all four folds.
positive_cv <- aggregate(delta_RMSE ~ group, held_out,
                         function(z) all(z > 0))
names(positive_cv)[2] <- "all_four_folds_improve"
eligible <- merge(count, positive_cv, by = "group", sort = FALSE)
eligible <- eligible[eligible$n_SNPs >= 25 & eligible$all_four_folds_improve, ]
eligible <- eligible[order(eligible$n_SNPs, decreasing = TRUE), ]
stopifnot(nrow(eligible) >= 2L)
main_groups <- eligible$group[1:2]
write.csv(data.frame(group = main_groups,
  criterion = "two largest groups with >=25 SNPs and positive climate-added held-out RMSE improvement in all four folds"),
  file.path(base, "39_v4_main_map_group_selection.csv"), row.names = FALSE)

# 7. Plot signed associations, frequencies and geographic scores ----

gmt <- readLines(file.path(base, "south_america_borders.gmt"), warn = FALSE)
coordinate <- grepl("^[0-9-]+[.]?[0-9]*[[:space:]]+[-0-9]", gmt)
border <- read.table(text = gmt[coordinate],
                     col.names = c("longitude", "latitude"))
border$segment <- cumsum(gmt == ">")[coordinate]
border <- border[is.finite(border$longitude) & is.finite(border$latitude), ]
gf <- read.csv(file.path(project,
  "results/nonlinear_present_landscapes/GF_controlled_mean_cumulative_curves.csv"))
gf <- gf[gf$predictor %in% c("bio7", "bio2", "bio4"), ]
gf$predictor <- toupper(gf$predictor)
stopifnot(all(is.finite(gf$mean_cumimp)))

deep <- "#32134F"; muted <- "#9071AD"; orange <- "#F28E2B"
paper <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", colour = deep, size = 12),
        plot.tag = element_text(face = "bold", colour = deep, size = 12),
        legend.position = "bottom", legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 9))

p_a <- ggplot(heat, aes(predictor, group_label, fill = standardized_beta)) +
  geom_tile(colour = "white", linewidth = .35) +
  scale_fill_gradient2(low = deep, mid = "white", high = orange, midpoint = 0,
                       name = "Mean LFMM2 effect\n(ALT SD / climate SD)") +
  labs(tag = "a", title = "Signed climatic profiles",
       x = NULL, y = "Dominant predictor and sign · SNPs/RAD loci") + paper

frequency$sign <- ifelse(grepl("\\+$", frequency$group), "+", "-")
bins$sign <- ifelse(grepl("\\+$", bins$group), "+", "-")
frequency$climate_label <- factor(climate_names[frequency$predictor],
                                  levels = climate_names)
bins$climate_label <- factor(climate_names[bins$predictor],
                             levels = climate_names)
p_b <- ggplot(frequency, aes(climate_original, mean_ALT_frequency,
                             colour = sign)) +
  geom_point(alpha = .15, size = .65) +
  geom_line(data = bins, aes(mean_climate, mean_ALT_frequency,
                             group = sign), linewidth = .8) +
  geom_point(data = bins, aes(mean_climate, mean_ALT_frequency), size = 1.65) +
  facet_wrap(~ climate_label, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = c("-" = deep, "+" = orange),
                      name = "LFMM2 group") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(tag = "b", title = "Observed climatic frequencies",
       x = "Observed climate (facet units)", y = "Mean ALT frequency") + paper

main_map <- map_data[map_data$group %in% main_groups, ]
main_map$group <- factor(main_map$group, levels = main_groups)
p_c <- ggplot() +
  geom_path(data = border, aes(longitude, latitude, group = segment),
            colour = "grey72", linewidth = .2) +
  geom_point(data = main_map,
             aes(longitude, latitude, colour = display_z, size = n_individuals),
             alpha = .95) +
  facet_wrap(~ group, nrow = 1) +
  coord_equal(xlim = c(-83, -34), ylim = c(-36, 12), expand = FALSE) +
  scale_colour_gradient2(low = deep, mid = "white", high = orange,
                         midpoint = 0, limits = c(-display_limit, display_limit),
                         oob = scales::squish) +
  scale_size_continuous(range = c(1.2, 3.5)) +
  labs(tag = "c", title = "Observed group geography",
       x = "Longitude", y = "Latitude",
       colour = "Within-group\nlocality z", size = "Individuals") + paper

p_d <- ggplot(gf, aes(ecdf, mean_cumimp, colour = predictor)) +
  geom_line(linewidth = .9) +
  scale_colour_manual(values = c(BIO2 = deep, BIO4 = muted, BIO7 = orange)) +
  labs(tag = "d", title = "Full-panel GF turnover",
       x = "Observed climatic gradient (ECDF)",
       y = "Mean cumulative split importance", colour = "Predictor") + paper

# 8. Export main and supplementary figures ----

panels <- list(A = p_a, B = p_b, C = p_c, D = p_d)
for (i in seq_along(panels)) {
  ggsave(file.path(base, paste0("Figure_candidate_directionality_v4_panel_",
                                names(panels)[i], ".pdf")),
         panels[[i]], width = 7.8, height = 5.3, bg = "white")
  ggsave(file.path(base, paste0("Figure_candidate_directionality_v4_panel_",
                                names(panels)[i], ".png")),
         panels[[i]], width = 7.8, height = 5.3, dpi = 400, bg = "white")
}
for (device in c("pdf", "png")) {
  path <- file.path(base, paste0("Figure_candidate_directionality_v4.", device))
  if (device == "pdf") pdf(path, width = 15.6, height = 10.6)
  if (device == "png") png(path, width = 15.6, height = 10.6,
                            units = "in", res = 400, bg = "white")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
  for (i in seq_along(panels)) {
    print(panels[[i]], vp = grid::viewport(layout.pos.row = (i + 1) %/% 2,
                                          layout.pos.col = (i - 1) %% 2 + 1))
  }
  grid::popViewport()
  dev.off()
}

p_frequency <- ggplot(frequency, aes(climate_original, mean_ALT_frequency,
                                     colour = sign)) +
  geom_point(aes(size = n_individuals), alpha = .32) +
  geom_line(data = bins, aes(mean_climate, mean_ALT_frequency,
                             group = sign), linewidth = .9) +
  geom_point(data = bins, aes(mean_climate, mean_ALT_frequency), size = 2) +
  facet_wrap(~ climate_label, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = c("-" = deep, "+" = orange),
                      name = "LFMM2 group") +
  scale_size_continuous(range = c(.6, 2.6), name = "Individuals") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "All eight observed frequency profiles",
       x = "Observed climate (facet units)", y = "Mean ALT frequency") + paper
ggsave(file.path(base, "Figure_all_eight_frequency_profiles_v4.pdf"),
       p_frequency, width = 11, height = 8, bg = "white")
ggsave(file.path(base, "Figure_all_eight_frequency_profiles_v4.png"),
       p_frequency, width = 11, height = 8, dpi = 400, bg = "white")

map_data$group <- factor(map_data$group, levels = count$group)
p_maps <- ggplot() +
  geom_path(data = border, aes(longitude, latitude, group = segment),
            colour = "grey75", linewidth = .17) +
  geom_point(data = map_data,
             aes(longitude, latitude, colour = display_z, size = n_individuals),
             alpha = .94) +
  facet_wrap(~ group, ncol = 2) +
  coord_equal(xlim = c(-83, -34), ylim = c(-36, 12), expand = FALSE) +
  scale_colour_gradient2(low = deep, mid = "white", high = orange,
                         midpoint = 0, limits = c(-display_limit, display_limit),
                         oob = scales::squish) +
  scale_size_continuous(range = c(.9, 2.8)) +
  labs(title = "Observed geography of all eight groups",
       x = "Longitude", y = "Latitude",
       colour = "Within-group\nlocality z", size = "Individuals") + paper
ggsave(file.path(base, "Figure_all_eight_group_maps_v4.pdf"),
       p_maps, width = 11, height = 18, bg = "white")
ggsave(file.path(base, "Figure_all_eight_group_maps_v4.png"),
       p_maps, width = 11, height = 18, dpi = 400, bg = "white")

held_out$group <- factor(held_out$group, levels = count$group)
p_diagnostic <- ggplot(held_out,
  aes(delta_RMSE, group, colour = training_direction_matches_LFMM2)) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
  geom_point(aes(shape = factor(held_out_fold)), size = 2.3,
             position = position_dodge(width = .45)) +
  scale_colour_manual(values = c("TRUE" = deep, "FALSE" = orange),
                      name = "LFMM2 sign retained") +
  labs(title = "Spatial holdout: climate contribution",
       x = "Held-out RMSE reduction from adding climate",
       y = "LFMM2 group", shape = "Fold") + paper
ggsave(file.path(base, "Figure_group_spatial_diagnostics_v4.pdf"),
       p_diagnostic, width = 8, height = 5.6, bg = "white")
ggsave(file.path(base, "Figure_group_spatial_diagnostics_v4.png"),
       p_diagnostic, width = 8, height = 5.6, dpi = 400, bg = "white")

rho <- cor(candidate_scores$candidate_PC1, PCA$PC1, method = "spearman")
stopifnot(is.finite(rho), abs(rho - 0.907428865188619) < 1e-8)
write.csv(data.frame(n_individuals = nrow(PCA), spearman_rho = rho),
          file.path(base, "40_v4_candidate_vs_genome_PC1_check.csv"),
          row.names = FALSE)
