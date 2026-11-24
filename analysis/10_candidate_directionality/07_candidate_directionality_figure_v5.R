###############################################################################
#### FIGURE 2: CANDIDATE DIRECTIONALITY AND GEOGRAPHY ####
###############################################################################

#### 1. VALIDATED FIGURE INPUTS ####

library(ggplot2)

base <- Sys.getenv("AMEIVA_DEVELOPMENT_DIR", unset = getwd())
project <- Sys.getenv("AMEIVA_PROJECT_DIR")
stopifnot(dir.exists(base), nzchar(project), dir.exists(project))

count <- read.csv(file.path(base, "30_v4_group_counts.csv"))
heat <- read.csv(file.path(base, "31_v4_group_mean_signed_effects.csv"))
frequency <- read.csv(file.path(base, "32_v4_all_group_locality_ALT_frequency.csv"))
bins <- read.csv(file.path(base, "33_v4_observed_frequency_by_climate_decile.csv"))
site_scores <- read.csv(file.path(base, "18_LFMM2_group_locality_scores_v2.csv"),
                        check.names = FALSE)
map_scale <- read.csv(file.path(base, "38_v4_map_display_scale.csv"))
map_selection <- read.csv(file.path(base, "39_v4_main_map_group_selection.csv"))
gf <- read.csv(file.path(project,
  "results/nonlinear_present_landscapes/GF_controlled_mean_cumulative_curves.csv"))

stopifnot(nrow(count) == 8L, nrow(heat) == 32L,
          nrow(frequency) == 1480L, nrow(bins) == 80L,
          nrow(site_scores) == 185L,
          identical(map_selection$group, c("BIO14 +", "BIO13 -")),
          all(is.finite(heat$standardized_beta)),
          all(is.finite(frequency$mean_ALT_frequency)),
          all(is.finite(bins$mean_ALT_frequency)))

#### 2. DISPLAY DATA; NO MODEL FITTING ####

heat$group_label <- factor(heat$group_label,
  levels = paste0(count$group, "  ", count$n_SNPs, "/", count$n_RAD_loci))
heat$predictor <- factor(heat$predictor,
  levels = c("BIO2", "BIO3", "BIO13", "BIO14"))
effect_limit <- max(abs(heat$standardized_beta))

climate_names <- c(BIO2 = "BIO2  Diurnal temperature range (°C)",
                   BIO3 = "BIO3  Isothermality (ratio)",
                   BIO13 = "BIO13  Wettest-month precipitation (mm)",
                   BIO14 = "BIO14  Driest-month precipitation (mm)")
frequency$sign <- ifelse(grepl("\\+$", frequency$group), "+", "−")
bins$sign <- ifelse(grepl("\\+$", bins$group), "+", "−")
frequency$climate_label <- factor(climate_names[frequency$predictor],
                                  levels = climate_names)
bins$climate_label <- factor(climate_names[bins$predictor],
                             levels = climate_names)

map_data <- do.call(rbind, lapply(map_selection$group, function(g)
  data.frame(site_scores[, c("longitude", "latitude", "n_individuals")],
             group = g, raw_locality_score = site_scores[[g]])))
map_data$display_z <- ave(map_data$raw_locality_score, map_data$group,
                          FUN = function(z) as.numeric(scale(z)))
display_limit <- map_scale$symmetric_display_limit[1]
stopifnot(nrow(map_data) == 370L, is.finite(display_limit),
          display_limit > 0,
          isTRUE(all.equal(display_limit,
            max(abs(quantile(do.call(c, lapply(count$group, function(g)
              as.numeric(scale(site_scores[[g]])))), c(.02, .98)))))))
map_data$group <- factor(map_data$group, levels = map_selection$group)

gmt <- readLines(file.path(base, "south_america_borders.gmt"), warn = FALSE)
coordinate <- grepl("^[0-9-]+[.]?[0-9]*[[:space:]]+[-0-9]", gmt)
border <- read.table(text = gmt[coordinate],
                     col.names = c("longitude", "latitude"))
border$segment <- cumsum(gmt == ">")[coordinate]
border <- border[is.finite(border$longitude) & is.finite(border$latitude), ]

gf <- gf[tolower(gf$predictor) %in% c("bio7", "bio2", "bio4"), ]
gf$predictor <- factor(toupper(gf$predictor),
                       levels = c("BIO7", "BIO2", "BIO4"))
stopifnot(nrow(gf) > 0L, all(is.finite(gf$ecdf)),
          all(is.finite(gf$mean_cumimp)),
          setequal(as.character(unique(gf$predictor)), c("BIO7", "BIO2", "BIO4")))

#### 3. MANUSCRIPT COLOURS AND PRINT-SIZE TYPOGRAPHY ####

deep <- "#32134F"
muted <- "#9071AD"
lavender <- "#B7A5CA"
orange <- "#F28E2B"
paper <- theme_classic(base_size = 9.5, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", colour = deep, size = 10.2,
                                  margin = margin(b = 3)),
        plot.tag = element_text(face = "bold", colour = deep, size = 10.5),
        plot.tag.position = c(.012, .96),
        plot.margin = margin(3, 5, 3, 5),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8.2, colour = "#30283A"),
        legend.position = "bottom",
        legend.title = element_text(size = 8.2),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "#F5F0F7", colour = NA),
        strip.text = element_text(size = 8.5, face = "bold", colour = deep),
        panel.spacing = grid::unit(3, "mm"))

#### 4. FOUR PRINT-READY PANELS ####

p_a <- ggplot(heat, aes(predictor, group_label, fill = standardized_beta)) +
  geom_tile(colour = "white", linewidth = .35) +
  scale_fill_gradient2(low = deep, mid = "white", high = orange,
    midpoint = 0, limits = c(-effect_limit, effect_limit),
    breaks = c(-.4, -.2, 0, .2, .4),
    name = "Mean LFMM2 effect\n(ALT SD / climate SD)") +
  guides(fill = guide_colorbar(title.position = "top", label.position = "bottom",
    barwidth = grid::unit(1.85, "in"), barheight = grid::unit(.13, "in"))) +
  labs(title = "a  Signed LFMM2 profiles",
       x = NULL, y = "Group · SNPs/RAD loci") + paper +
  theme(axis.text.y = element_text(size = 8.7),
        axis.text.x = element_text(size = 8.8))

p_b <- ggplot(frequency, aes(climate_original, mean_ALT_frequency,
                             colour = sign)) +
  geom_point(alpha = .27, size = .82, show.legend = FALSE) +
  geom_line(data = bins, aes(mean_climate, mean_ALT_frequency, group = sign),
            linewidth = .9) +
  geom_point(data = bins, aes(mean_climate, mean_ALT_frequency),
             size = 1.25, show.legend = FALSE) +
  facet_wrap(~ climate_label, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = c("−" = deep, "+" = orange),
    breaks = c("−", "+"), labels = c("Negative", "Positive"),
    name = "LFMM2 ALT association") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, .5, 1)) +
  guides(colour = guide_legend(nrow = 1,
    override.aes = list(linewidth = 1.2, alpha = 1))) +
  labs(title = "b  Observed ALT frequencies",
       x = NULL, y = "Mean ALT frequency") + paper +
  theme(strip.text = element_text(size = 8.8, face = "bold", colour = deep))

p_c <- ggplot() +
  geom_path(data = border, aes(longitude, latitude, group = segment),
            colour = "#B6B0B9", linewidth = .19) +
  geom_point(data = map_data,
    aes(longitude, latitude, fill = display_z, size = n_individuals),
    shape = 21, colour = "#655874", stroke = .13, alpha = .98) +
  facet_wrap(~ group, nrow = 1) +
  coord_equal(xlim = c(-83, -34), ylim = c(-36, 12), expand = FALSE) +
  scale_fill_gradient2(low = deep, mid = "white", high = orange,
    midpoint = 0, limits = c(-display_limit, display_limit),
    oob = scales::squish, breaks = c(-2, -1, 0, 1, 2),
    name = "Standardized multilocus score") +
  scale_size_continuous(range = c(1.2, 3.2), breaks = 1:3,
                        name = "Individuals / locality") +
  guides(fill = guide_colorbar(order = 1, title.position = "top",
           barwidth = grid::unit(1.65, "in"),
           barheight = grid::unit(.13, "in")),
         size = guide_legend(order = 2, title.position = "top", nrow = 1)) +
  labs(title = "c  Sampled locality scores",
       x = "Longitude", y = "Latitude") + paper +
  theme(legend.box = "horizontal", legend.box.just = "center",
        legend.spacing.x = grid::unit(3, "mm"))

p_d <- ggplot(gf, aes(ecdf, mean_cumimp,
                       colour = predictor, linetype = predictor)) +
  geom_line(linewidth = .9) +
  scale_colour_manual(values = c(BIO7 = deep, BIO2 = muted,
                                  BIO4 = lavender), name = "Predictor") +
  scale_linetype_manual(values = c(BIO7 = "solid", BIO2 = "longdash",
                                   BIO4 = "dotted"), name = "Predictor") +
  scale_x_continuous(breaks = seq(0, 1, .25)) +
  labs(title = "d  Full-panel GF turnover",
       x = "Observed climate ECDF",
       y = "Cumulative split importance") + paper

#### 5. FINAL TWO-COLUMN EXPORT ####

width_mm <- 176
height_mm <- 234
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    3, 2, widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(2.4, 2.8, 4.0), "null"))))
  print(p_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_d, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p_b, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  print(p_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2))
  grid::popViewport()
}

pdf_file <- file.path(base, "Figure_candidate_directionality_v5.pdf")
png_file <- file.path(base, "Figure_candidate_directionality_v5.png")
grDevices::cairo_pdf(pdf_file, width = width_in, height = height_in,
                     family = "sans", bg = "white")
draw_figure()
dev.off()
grDevices::png(png_file, width = width_in, height = height_in,
               units = "in", res = 400, type = "cairo", bg = "white")
draw_figure()
dev.off()

message("Figure 2 exported at ", width_mm, " × ", height_mm,
        " mm; all data and model outputs unchanged.")
