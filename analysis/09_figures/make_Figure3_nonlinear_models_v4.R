
# 1. Read fitted GF summaries, locus curves and RDAforest importance ----

library(ggplot2)
library(patchwork)

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"
base_results <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation/results/RDA-2027-14134"
gf_dir <- file.path(base_results, "gradient_forest_all_14134")
followup <- file.path(project, "results", "GF_locus_followup")
main_figures <- file.path(project, "results", "figures", "RDA-2027-14134")
manuscript_figures <- file.path(project, "manuscript_figures_final")
environment_file <- "/home/nadaline/Documents/Phd/RADseq/Ameiva_Variables_rda.csv"

gf_importance <- read.csv(file.path(gf_dir, "controlled_PC1_PC3_lon_lat",
                                    "GF_controlled_predictor_importance_by_block.csv"))
locus <- read.csv(file.path(followup, "01_GF_locus_summary.csv"))
turnover_curves <- read.csv(file.path(followup, "04_selected_locus_extracted_turnover_curves.csv"))
rdaforest <- read.csv(file.path(base_results, "12_RDAforest_importance.csv"))
variables <- read.csv(environment_file, check.names = FALSE, stringsAsFactors = FALSE)

predictors <- c("BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6", "BIO7", "BIO12", "BIO13", "BIO15")
high <- locus[locus$GF_high, ]
stopifnot(nrow(locus) == 14134L, nrow(high) == 707L)

# 2. Define figure colours and theme ----

deep_purple <- "#32134F"
purple <- "#5B2A86"
muted_purple <- "#9071AD"
lavender <- "#D8C9E5"
pale_warm <- "#F5E3DA"
orange <- "#F28E2B"

panel_theme <- theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", colour = deep_purple, size = 10.5, margin = margin(b = 5)),
    plot.tag = element_text(face = "bold", colour = deep_purple, size = 10),
    axis.title = element_text(colour = "grey20", size = 9),
    axis.text = element_text(colour = "grey20", size = 8.5),
    axis.line = element_line(linewidth = .3),
    axis.ticks = element_line(linewidth = .3),
    legend.title = element_text(face = "bold", colour = deep_purple, size = 8.5),
    legend.text = element_text(size = 8)
  )

# 3. Summarise predictor importance across Gradient Forest batches ----

gf_summary <- aggregate(importance ~ predictor, gf_importance,
                        function(x) c(mean = mean(x), sd = sd(x)))
gf_summary <- data.frame(
  predictor = toupper(gf_summary$predictor),
  mean_importance = gf_summary$importance[, "mean"],
  sd_importance = gf_summary$importance[, "sd"]
)
gf_summary <- gf_summary[order(gf_summary$mean_importance, decreasing = TRUE), ]

p_a <- ggplot(gf_summary, aes(reorder(predictor, mean_importance), mean_importance,
                              fill = mean_importance)) +
  geom_col(width = .68) +
  coord_flip() +
  scale_fill_gradientn(colours = c(pale_warm, lavender, muted_purple, purple, deep_purple), guide = "none") +
  labs(tag = "a", title = "Gradient Forest predictor importance",
       x = NULL, y = "Mean GF importance") +
  panel_theme

# 4. Count high-performing loci by predictor and turnover position ----

high$turnover_position <- ifelse(high$peak_split_environment_quantile < .33, "Low",
                                 ifelse(high$peak_split_environment_quantile <= .67,
                                        "Intermediate", "High"))
high$turnover_position <- factor(high$turnover_position, levels = c("Low", "Intermediate", "High"))

mosaic <- as.data.frame(table(high$top_predictor, high$turnover_position), stringsAsFactors = FALSE)
names(mosaic) <- c("predictor", "gradient_position", "n_loci")
mosaic$total <- ave(mosaic$n_loci, mosaic$predictor, FUN = sum)
mosaic$predictor <- as.character(mosaic$predictor)
# Use the same predictor order for the aggregate and locus-specific summaries.
predictor_order <- as.character(gf_summary$predictor)
mosaic$predictor <- factor(mosaic$predictor, levels = rev(predictor_order))
bar_ends <- mosaic[match(predictor_order, mosaic$predictor), c("predictor", "total")]
bar_ends$predictor <- factor(bar_ends$predictor, levels = rev(predictor_order))

p_b <- ggplot(mosaic, aes(predictor, n_loci, fill = gradient_position)) +
  geom_col(width = .70, colour = "white", linewidth = .16) +
  geom_text(data = bar_ends, aes(x = predictor, y = total, label = total), inherit.aes = FALSE,
            hjust = -.25, colour = deep_purple, size = 2.9, fontface = "bold") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Low" = deep_purple, "Intermediate" = lavender, "High" = orange),
                    name = "Position along observed\nclimatic gradient",
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, .10))) +
  labs(tag = "b", title = "Locus-specific climatic mosaic", x = NULL, y = "High-GF loci") +
  panel_theme +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(t = 3))

# 5. Plot four GF-only locus turnover curves ----

selected_ids <- c("SNP03109", "SNP06251", "SNP07293", "SNP12662")
selected_meta <- high[match(selected_ids, high$SNP), ]
stopifnot(identical(selected_meta$SNP, selected_ids), all(selected_meta$GF_only),
          all(selected_meta$checkerboard_all_folds))

environment_columns <- paste0("X.CHELSA_", tolower(predictors), "_1981.2010_SA")
environment_raw <- variables[, environment_columns, drop = FALSE]
names(environment_raw) <- tolower(predictors)
predictor_labels <- c(
  bio1 = "BIO1: annual mean temperature", bio2 = "BIO2: mean diurnal range",
  bio3 = "BIO3: isothermality", bio4 = "BIO4: temperature seasonality",
  bio5 = "BIO5: maximum temperature of warmest month",
  bio6 = "BIO6: minimum temperature of coldest month",
  bio7 = "BIO7: annual temperature range", bio12 = "BIO12: annual precipitation",
  bio13 = "BIO13: precipitation of wettest month", bio15 = "BIO15: precipitation seasonality"
)

curve_data <- do.call(rbind, lapply(seq_len(nrow(selected_meta)), function(i) {
  snp <- selected_meta$SNP[i]
  predictor <- tolower(selected_meta$top_predictor[i])
  d <- turnover_curves[turnover_curves$SNP == snp & turnover_curves$var == predictor, ]
  stopifnot(nrow(d) > 0)
  d$x_raw <- as.numeric(quantile(environment_raw[[predictor]], d$split, type = 8, names = FALSE))
  d$SNP <- snp
  d$predictor_label <- unname(predictor_labels[predictor])
  d
}))
curve_data$SNP <- factor(curve_data$SNP, levels = selected_ids)
curve_data$panel_label <- paste0(curve_data$SNP, "\n", curve_data$predictor_label)
curve_data$panel_label <- factor(curve_data$panel_label,
                                 levels = unique(curve_data$panel_label)[match(selected_ids, unique(curve_data$SNP))])

# Mark the split with the largest improvement in each displayed curve.
curve_peaks <- do.call(rbind, lapply(split(curve_data, curve_data$SNP), function(d) d[which.max(d$improve), ]))

p_c <- ggplot(curve_data, aes(x_raw, cumulative_fraction)) +
  geom_line(colour = deep_purple, linewidth = .72) +
  geom_vline(data = curve_peaks, aes(xintercept = x_raw), inherit.aes = FALSE,
             colour = orange, linewidth = .55, linetype = "dashed") +
  facet_wrap(~ panel_label, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, .02))) +
  labs(tag = "c", title = "Spatially represented locus-specific turnover",
       x = "Observed climatic predictor (original units)", y = "Cumulative GF split importance") +
  panel_theme +
  theme(strip.background = element_rect(fill = "white", colour = "grey78", linewidth = .25),
        strip.text = element_text(face = "bold", size = 8, colour = deep_purple),
        panel.spacing = grid::unit(8, "pt"))

# 6. Summarise RDAforest predictor importance ----

rdaforest$variable <- factor(toupper(rdaforest$variable), levels = c("BIO3", "BIO13", "BIO2"))
rdaforest_mean <- aggregate(importance ~ variable, rdaforest, mean)
p_d <- ggplot(rdaforest, aes(variable, importance)) +
  geom_boxplot(width = .56, fill = lavender, colour = deep_purple, linewidth = .4, outlier.shape = NA) +
  geom_point(position = position_jitter(width = .08, height = 0), size = .8, alpha = .48, colour = muted_purple) +
  geom_point(data = rdaforest_mean, aes(variable, importance), inherit.aes = FALSE,
             colour = orange, size = 1.75) +
  coord_flip() +
  labs(tag = "d", title = "RDAforest predictor contribution", x = NULL,
       y = "Proportion of total genomic variation explained") +
  panel_theme

# 7. Assemble and export the nonlinear model figure ----

figure_3_v4 <- (p_a + p_d) / p_b / p_c + plot_layout(heights = c(1.25, 1.65, 2.20))

for (extension in c("pdf", "png")) {
  file <- file.path(followup, paste0("Figure_3_nonlinear_models_v4.", extension))
  if (extension == "pdf") {
    ggsave(file, figure_3_v4, width = 12, height = 13.2, device = cairo_pdf)
  } else {
    ggsave(file, figure_3_v4, width = 12, height = 13.2, dpi = 600, bg = "white")
  }
}
ggsave(file.path(main_figures, "Main_Figure_3_nonlinear_models.png"), figure_3_v4,
       width = 12, height = 13.2, dpi = 600, bg = "white")
ggsave(file.path(main_figures, "Main_Figure_3_nonlinear_models.tiff"), figure_3_v4,
       width = 12, height = 13.2, dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(manuscript_figures, "Figure_3_nonlinear_models.png"), figure_3_v4,
       width = 12, height = 13.2, dpi = 600, bg = "white")
ggsave(file.path(manuscript_figures, "Figure_3_nonlinear_models.tiff"), figure_3_v4,
       width = 12, height = 13.2, dpi = 600, compression = "lzw", bg = "white")

# 8. Export supplementary locus and turnover figures ----

# The supplementary figures show all 707 loci and the individual-level diagnostics
# for the four illustrative SNPs.
file.copy(file.path(followup, "Figure_3_nonlinear_models_v2.png"),
          file.path(followup, "Supplementary_Figure_GF_707_locus_detail.png"), overwrite = TRUE)
file.copy(file.path(followup, "Figure_GF_four_spatially_represented_loci.png"),
          file.path(followup, "Supplementary_Figure_GF_four_locus_diagnostics.png"), overwrite = TRUE)

caption <- paste(
  "Figure 3. Nonlinear climate-associated genomic mosaic in Ameiva ameiva.",
  "(a) Mean controlled Gradient Forest (GF) importance across 32 computational subsets.",
  "(b) Among the 707 high-GF loci (top 5% of controlled OOB R2), the climatic predictor with the highest locus-specific importance is distributed across low, intermediate and high portions of its observed ECDF-scaled gradient.",
  "These descriptive positions are neither effect directions nor adaptive thresholds.",
  "(c) Directly extracted cumulative GF turnover curves for four GF-only loci retaining residual variation in all four checkerboard folds; axes are shown in original climatic units and dashed orange lines mark the largest stored split improvement within each displayed top-predictor curve.",
  "(d) Jackknife distributions of predictor contributions retained in the main RDAforest model.",
  "GF quantities are descriptive model attributions; the four loci illustrate, rather than exhaust, locus-level response shapes."
)
writeLines(caption, file.path(followup, "Figure_3_caption_draft_v4.txt"))

message("Figure 3 v4 complete | synthesis only | no models or values recalculated")
