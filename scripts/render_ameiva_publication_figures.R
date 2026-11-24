# 1. Load packages and define shared export functions ----

library(ggplot2)
library(patchwork)
library(maps)

project <- normalizePath(getwd())
root <- Sys.getenv("AMEIVA_RESULTS_DIR", "results/RDA-2027-14134")
out <- Sys.getenv("AMEIVA_PUBLICATION_DIR", "results/publication_figures_2026-11-19")
stopifnot(dir.exists(root), dir.exists(out))
purple <- "#32134F"
orange <- "#F28E2B"
cream <- "#FFF7EF"

print_theme <- theme(axis.text = element_text(size = 9),
  axis.title = element_text(size = 10), legend.title = element_text(size = 9),
  legend.text = element_text(size = 9), strip.text = element_text(size = 9),
  plot.title = element_text(size = 11, face = "bold", colour = purple),
  plot.tag = element_text(size = 13, face = "bold", colour = purple))

export_plate <- function(plot, stem, height, width = 178) {
  paths <- file.path(out, paste0(stem, c(".pdf", ".png", ".tiff")))
  stopifnot(!any(file.exists(paths)))
  ggsave(paths[1], plot, width = width, height = height, units = "mm", bg = "white", device = cairo_pdf)
  ggsave(paths[2], plot, width = width, height = height, units = "mm", bg = "white", dpi = 600, type = "cairo")
  ggsave(paths[3], plot, width = width, height = height, units = "mm", bg = "white", dpi = 600,
         type = "cairo", compression = "lzw")
}

# Read the plot-building expressions up to the designated export command.
read_plot_code <- function(path, stop_at, skip = character()) {
  lines <- readLines(path, warn = FALSE)
  j <- which(grepl(stop_at, lines, fixed = TRUE))[1]
  stopifnot(is.finite(j))
  lines <- lines[seq_len(j - 1)]
  if (length(skip)) lines <- lines[!Reduce("|", lapply(skip, grepl, x = lines, fixed = TRUE))]
  eval(parse(text = lines), envir = new.env(parent = globalenv()))
}

tasks <- commandArgs(TRUE)
if (!length(tasks)) tasks <- c("1", "2", "3", "4", "5", "offset_rank", "S1", "demography", "S7", "directionality", "diagnostics", "frequency", "offset_supplement")

read_xyz <- function(path, band = 1L) {
  con <- pipe(paste("gdal_translate -q -of XYZ -b", band, shQuote(path), "/vsistdout/"), open = "r")
  on.exit(close(con))
  x <- read.table(con, col.names = c("x", "y", "score"), na.strings = c("NA", "nan"))
  dx <- median(diff(sort(unique(x$x))))
  dy <- median(diff(sort(unique(x$y))))
  x$x <- round((x$x - x$x[1])/dx)*dx + x$x[1]
  x$y <- round((x$y - x$y[1])/dy)*dy + x$y[1]
  x
}

# 2. Render population structure, linear associations and nonlinear synthesis ----

if ("1" %in% tasks) {
  src <- readLines("scripts/restyle_main_figures_1_2.R")
  src <- src[seq_len(which(startsWith(src, "save_plate(figure_1,")) - 1L)]
  src <- src[!src %in% c("library(vegan)", "library(terra)")]
  src <- sub("Genomic structure across Ameiva ameiva", "Genome-wide PCA", src, fixed = TRUE)
  src <- sub('geom_col(width = 1.01, colour = NA)', 'geom_col(width = 1, colour = NA)', src, fixed = TRUE)
  src <- sub('scale_y_continuous(expand = c(0, 0), limits = c(0, 1))',
    'scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(0, 1))', src, fixed = TRUE)
  f1 <- eval(parse(text = src), new.env(parent = globalenv()))
  export_plate(f1 & print_theme, "Figure_1_population_structure", 170)
}

if ("2" %in% tasks) {
  if (nzchar(Sys.getenv("AMEIVA_VEGAN_LIBRARY")))
    .libPaths(c(.libPaths(), Sys.getenv("AMEIVA_VEGAN_LIBRARY")))
  library(vegan)
  lines <- readLines("scripts/restyle_main_figures_1_2.R")
  first <- which(startsWith(lines, "# Existing present-day RDA axes"))
  plot_a <- which(startsWith(lines, "p2_landscape <-"))
  df_b <- which(startsWith(lines, "rda2_landscape_df <-"))
  plot_b <- which(startsWith(lines, "p2_landscape_rda2 <-"))
  last <- which(startsWith(lines, "save_plate(figure_2,"))
  e <- new.env()
  init <- lines[seq_len(first - 1L)]
  init <- init[!grepl('library\\(terra\\)|^save_plate\\(figure_1|^file.copy', init)]
  eval(parse(text = init), e)
  raw <- lapply(1:2, function(i) read_xyz(file.path(root, "spatial_offset/18_RDA_present_axes.tif"), i))
  for (i in 1:2) {
    d <- raw[[i]]
    lim <- max(abs(quantile(d$score[is.finite(d$score)], c(.02, .98)))) * c(-1, 1)
    # Average raster cells in 4 × 4 blocks for the display.
    dx <- median(diff(sort(unique(d$x))))
    dy <- median(diff(sort(unique(d$y))))
    west <- min(d$x); north <- max(d$y)
    ix <- round((d$x-west)/dx) %/% 4L
    iy <- round((north-d$y)/dy) %/% 4L
    ok <- is.finite(d$score)
    df <- aggregate(score ~ ix + iy, data.frame(ix, iy, score = d$score)[ok, ], mean)
    df <- data.frame(x = west + (4*df$ix + 1.5)*dx,
                     y = north - (4*df$iy + 1.5)*dy, score = df$score)
    names(df)[3] <- paste0("RDA", i)
    assign(c("rda_landscape_df", "rda2_landscape_df")[i], df, e)
    assign(paste0("rda", i, "_display_limit"), lim, e)
  }
  eval(parse(text = c(lines[plot_a:(df_b-1L)], lines[plot_b:(last-1L)])), e)
  small_bar <- guides(fill = guide_colourbar(title.position = "top",
    barwidth = unit(33, "mm"), barheight = unit(2.6, "mm")))
  e$p2_rda$layers[[5]]$data$environment <- toupper(e$p2_rda$layers[[5]]$data$environment)
  p <- ((e$p2_rda + labs(title = "Structure-controlled pRDA") |
        e$p2_pca + labs(title = "Linear candidate clusters")) /
        (e$p2_map + labs(title = "Candidate clusters") |
         e$p2_landscape + labs(title = "RDA1 landscape") + small_bar |
         e$p2_landscape_rda2 + labs(title = "RDA2 landscape") + small_bar)) +
    plot_layout(heights = c(1, .95)) + plot_annotation(tag_levels = "a")
  export_plate(p & print_theme, "Figure_2_linear_climate_associations", 173)
}

if ("3" %in% tasks) {
  f3 <- read_plot_code("scripts/make_Figure3_nonlinear_models_v4.R", 'for (extension')
  export_plate(f3 & print_theme, "Figure_3_nonlinear_mosaic", 226)
}

# 3. Render present-day genomic landscapes ----

if ("4" %in% tasks) {
  landscape <- "results/nonlinear_present_landscapes"
  xyz <- function(name) {
    x <- read_xyz(file.path(landscape, name))
    meta <- system2("gdalinfo", shQuote(file.path(landscape, name)), stdout = TRUE)
    nd <- grep("NoData Value=", meta, value = TRUE)
    if (length(nd)) {
      nodata <- suppressWarnings(as.numeric(sub(".*NoData Value=", "", nd[1])))
      if (is.finite(nodata)) x$score[x$score == nodata] <- NA_real_
    }
    x[is.finite(x$score), ]
  }
  world <- map_data("world")
  samples <- read.csv(Sys.getenv("AMEIVA_ENV_FILE", "../Ameiva_Variables_rda.csv"))
  direction <- read.csv(file.path(landscape, "Figure_4_PC1_climate_direction.csv"))
  # Draw the saved landscape scores and climate-correlation insets.
  lines <- readLines("scripts/finalize_Figure4_present_nonlinear_landscapes.R")
  first <- which(startsWith(lines, "  ggplot(d, aes(x, y, fill = score))"))
  last <- which(startsWith(lines, "# RDAforest has")) - 2L
  map_lines <- sub("geom_raster() +", "geom_tile(width = dx, height = dy) +", lines[first:last], fixed = TRUE)
  map_saved <- eval(parse(text = c(
    "function(d, native, title, legend) {",
    "lim <- max(abs(quantile(native$score, c(.02, .98)))) * c(-1, 1)",
    "dx <- min(diff(sort(unique(d$x)))); dy <- min(diff(sort(unique(d$y))))",
    map_lines)))
  inset_first <- which(startsWith(lines, "climate_inset <-"))
  inset_last <- which(startsWith(lines, "inset_gf <-")) - 1L
  eval(parse(text = lines[inset_first:inset_last]))
  gf <- xyz("GF_present_landscape.tif")
  rdf <- xyz("RDAforest_present_landscape.tif")
  bilinear <- xyz("RDAforest_present_landscape_bilinear_display_only.tif")
  p_gf <- map_saved(gf, gf, "GF adaptive landscape", "GF PC1") + labs(tag = "a")
  p_rdf <- map_saved(bilinear, rdf, "RDAforest-derived\nadaptive landscape", "RDAforest PC1") + labs(tag = "b")
  p_i <- climate_inset(head(subset(direction, Method == "Gradient Forest"), 4)) +
    labs(y = expression("Spearman " * rho))
  p_j <- climate_inset(subset(direction, Method == "RDAforest")) +
    labs(y = expression("Spearman " * rho))
  export_plate(((p_gf | p_rdf) / (p_i | p_j) + plot_layout(heights = c(3, 1.55))) &
    print_theme, "Figure_4_present_adaptive_landscapes", 158)
}

# 4. Compare mid- and late-century genomic offset ranks ----

if (any(c("5", "offset_rank", "s6_display", "s10_magnitudes") %in% tasks)) {
  x <- read.csv("results/Figure_5_genomic_offset_final/Figure_5_relative_offset_source.csv")
  mask <- read.csv(file.path(out, "gf_late_grid_mask.csv"))
  late <- subset(x, period == "2071-2100")
  mid <- subset(x, period == "2041-2070")
  stopifnot(isTRUE(all.equal(late[c("x", "y")], mask[c("x", "y")], check.attributes = FALSE)),
    identical(late$GF_valid, as.logical(mask$offset_finite)),
    all(late$GF_valid <= as.logical(mask$combined_supported)))
  late$novel <- as.logical(mask$future_novel)
  mid$novel <- FALSE
  stopifnot(!any(late$novel & late$GF_valid))
  rank01 <- function(v) (rank(v, ties.method = "average") - 1) / (length(v) - 1)
  eq <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12))
  for (d in list(mid, late)) {
    for (method in c("RDA", "GF", "RDAforest")) {
      valid <- d[[paste0(method, "_valid")]]
      ranks <- d[[paste0(method, "_rank")]]
      stopifnot(identical(valid, is.finite(d[[method]])),
        identical(valid, is.finite(ranks)),
        eq(ranks[valid], rank01(d[[method]][valid])))
      message(d$period[1], " ", method, " ranked: ", sum(valid),
              " of ", nrow(d), " cells")
    }
    shared <- d$RDA_valid & d$GF_valid & d$RDAforest_valid
    stopifnot(identical(shared, d$shared_valid), all(is.na(d$agreement[!shared])))
    linear <- rank01(d$RDA[shared])
    gf <- rank01(d$GF[shared])
    rdf <- rank01(d$RDAforest[shared])
    nonlinear <- (gf + rdf) / 2
    stopifnot(eq(linear, d$linear_shared_rank[shared]), eq(gf, d$GF_shared_rank[shared]),
      eq(rdf, d$RDAforest_shared_rank[shared]), eq(nonlinear, d$nonlinear_shared_rank[shared]))
    # Check agreement classes against the saved classifications.
    class_check <- rep("Intermediate", sum(shared))
    class_check[linear >= .67 & nonlinear >= .67] <- "High–high"
    class_check[linear <= .33 & nonlinear <= .33] <- "Low–low"
    class_check[linear - nonlinear > .25] <- "Linear higher"
    class_check[nonlinear - linear > .25] <- "Nonlinear higher"
    stopifnot(identical(class_check, d$agreement[shared]))
    message(d$period[1], " agreement verified: ", sum(shared), " shared cells")
  }
  relative_colours <- c("#FFF7EF", "#F8C58F", "#D99A8B", "#8F668D", "#5B2A86", "#32134F")
  world <- map_data("world")
  map_context <- geom_polygon(data = world, aes(long, lat, group = group),
    inherit.aes = FALSE, fill = NA, colour = "grey67", linewidth = .18)
  coordinate <- coord_quickmap(xlim = c(-85, -30), ylim = c(-35, 15), expand = FALSE)
  map_theme <- theme_minimal() + theme(panel.grid = element_blank(),
    legend.position = "bottom", legend.box = "vertical") + print_theme
  hx <- median(diff(sort(unique(late$x))))
  hy <- median(diff(sort(unique(late$y))))
  novelty_edges <- function(d) {
    # Trace the native cells of the binary domain mask.
    k <- function(x, y) paste(round(x / hx - min(d$x) / hx),
                               round(y / hy - min(d$y) / hy))
    n <- d[d$novel, c("x", "y")]
    keys <- k(n$x, n$y)
    do.call(rbind, lapply(list(c(-1,0), c(1,0), c(0,-1), c(0,1)), function(s) {
      q <- n[!k(n$x + s[1]*hx, n$y + s[2]*hy) %in% keys, , drop = FALSE]
      data.frame(x = q$x + (s[1]-abs(s[2]))*hx/2,
        xend = q$x + (s[1]+abs(s[2]))*hx/2,
        y = q$y + (s[2]-abs(s[1]))*hy/2,
        yend = q$y + (s[2]+abs(s[1]))*hy/2,
        domain = rep("Outside GF training climate", nrow(q)))
    }))
  }
  contour_layer <- function(d) {
    edges <- novelty_edges(d)
    list(geom_segment(data = edges, mapping = aes(x, y, xend = xend, yend = yend, colour = domain),
           inherit.aes = FALSE, linewidth = .50, linetype = "solid"),
         scale_colour_manual(values = c("Outside GF training climate" = "white"), name = NULL,
           guide = guide_legend(theme = theme(legend.key = element_rect(fill = purple, colour = NA)))))
  }
  offset_panel <- function(d, method, title, tag, relative, completed = FALSE) {
    variable <- if (relative) paste0(method, "_rank") else method
    scale_limits <- if (relative) c(0, 1) else range(x[[method]], finite = TRUE)
    colours <- if (relative) relative_colours else c(cream, "#F8C58F", "#D99A8B", "#8F668D", purple)
    use_display <- completed && method == "GF"
    if (use_display) variable <- "GF_display_rank"
    cells <- if (use_display) d else d[d[[paste0(method, "_valid")]], ]
    p <- ggplot(cells, aes(x, y, fill = .data[[variable]])) +
      geom_tile(width = hx * 1.01, height = hy * 1.01)
    if (use_display) {
      p <- p + contour_layer(d)
    } else if (method == "GF" && any(d$novel)) {
      p <- p + geom_tile(data = subset(d, novel), aes(x, y), inherit.aes = FALSE,
        width = hx * 1.01, height = hy * 1.01, fill = "#ECE8EF") +
        geom_segment(data = novelty_edges(d), aes(x, y, xend = xend, yend = yend, colour = domain),
          inherit.aes = FALSE, linewidth = .37) +
        scale_colour_manual(values = c("Outside GF training climate" = "#705185"), name = NULL)
    }
    p + map_context + coordinate +
      scale_fill_gradientn(colours = colours, limits = scale_limits,
        breaks = if (relative) c(0, .2, .4, .6, .8, 1) else waiver(),
        n.breaks = 4, name = if (relative) "Relative genomic offset" else "Genomic offset",
        guide = guide_colourbar(title.position = "top", barwidth = unit(34, "mm"))) +
      labs(title = title, tag = tag, x = "Longitude", y = "Latitude") + map_theme
  }
  agreement <- function(d, completed = FALSE) {
    stopifnot(all(is.na(d$agreement[!d$shared_valid])))
    cells <- if (completed) d else subset(d, shared_valid)
    variable <- if (completed) "display_agreement" else "agreement"
    p <- ggplot(cells, aes(x, y, fill = .data[[variable]])) +
      geom_tile(width = hx * 1.01, height = hy * 1.01)
    if (completed) p <- p + contour_layer(d)
    p + map_context + coordinate +
      scale_fill_manual(name = NULL,
        values = c("High–high" = purple, "Intermediate" = "#B6A2C8",
          "Low–low" = "#FFF1DF", "Linear higher" = "#73539B", "Nonlinear higher" = orange),
        guide = guide_legend(ncol = 2, byrow = TRUE)) +
      labs(title = "Spatial agreement", tag = "d", x = "Longitude", y = "Latitude") + map_theme
  }
  if ("s10_magnitudes" %in% tasks) {
    methods <- c("RDA", "GF", "RDAforest")
    titles <- c("Candidate-enriched RDA", "Gradient Forest", "RDAforest")
    s10_maps <- list()
    for (d in list(mid, late)) {
      for (j in seq_along(methods)) {
        m <- methods[j]
        valid <- d[[paste0(m, "_valid")]]
        limits <- range(x[[m]], finite = TRUE)
        has_novelty <- m == "GF" && any(d$novel)
        q <- ggplot(d[valid, ], aes(x, y, fill = .data[[m]])) +
          geom_tile(width = hx * 1.01, height = hy * 1.01)
        if (has_novelty) {
          # Mark cells without a supported GF estimate using a neutral background.
          q <- q + geom_tile(data = subset(d, novel), aes(x, y),
            inherit.aes = FALSE, width = hx * 1.01, height = hy * 1.01, fill = "#A394AC")
        }
        q <- q + map_context
        if (has_novelty) q <- q + contour_layer(d)
        q <- q + coordinate +
          scale_x_continuous(breaks = c(-80, -60, -40)) +
          scale_y_continuous(breaks = c(-30, -10, 10)) +
          scale_fill_gradientn(colours = relative_colours, limits = limits,
            name = "Original genomic offset",
            breaks = scales::breaks_extended(n = 4),
            guide = guide_colourbar(title.position = "top", barwidth = unit(35, "mm"),
                                    barheight = unit(3.4, "mm"))) +
          labs(title = if (m == "RDA") "RDA" else titles[j], subtitle = sub("-", "–", d$period[1]),
            tag = letters[length(s10_maps) + 1L], x = "Longitude", y = "Latitude") +
          map_theme + theme(legend.text = element_text(size = 10),
                           legend.title = element_text(size = 10),
                           legend.key.height = unit(3.3, "mm"))
        s10_maps[[length(s10_maps) + 1L]] <- q
      }
    }
    shared <- subset(x, shared_valid)
    stopifnot(identical(as.integer(table(shared$period)), c(3971L, 2529L)))
    histograms <- lapply(seq_along(methods), function(j) {
      m <- methods[j]
      limits <- range(x[[m]], finite = TRUE)
      breaks <- seq(limits[1], limits[2], length.out = 23)
      ggplot(shared, aes(.data[[m]], fill = period, colour = period)) +
        geom_histogram(aes(y = after_stat(density)), breaks = breaks,
                       position = "identity", alpha = .35, linewidth = .28) +
        coord_cartesian(xlim = limits) +
        scale_fill_manual(values = c("2041-2070" = purple, "2071-2100" = orange),
          labels = c("2041–2070", "2071–2100"), name = NULL) +
        scale_colour_manual(values = c("2041-2070" = purple, "2071-2100" = orange),
          labels = c("2041–2070", "2071–2100"), name = NULL) +
        labs(title = titles[j], x = "Original genomic offset", y = "Density",
          tag = letters[6L + j]) + theme_classic() + print_theme +
        theme(legend.position = "bottom", legend.text = element_text(size = 10))
    })
    maps <- wrap_plots(s10_maps, ncol = 3)
    distributions <- wrap_plots(histograms, ncol = 3, guides = "collect") &
      theme(legend.position = "bottom")
    export_plate(maps / distributions + plot_layout(heights = c(2.5, 1)),
                 "Figure_S10_original_offset_magnitudes_white_contour", 242, width = 190)
  }
  if ("s6_display" %in% tasks) {
    d <- late
    shared <- d$shared_valid
    filled <- !shared
    # Read the thin-plate GAM ranks used for the late-century display.
    stopifnot(all(is.finite(d$GF_display_rank)), all(d$GF_display_rank >= 0 & d$GF_display_rank <= 1),
      identical(d$GF_display_rank[d$GF_valid], d$GF_rank[d$GF_valid]))
    d$RDA_display_rank <- d$RDA_rank
    missing_rda <- !d$RDA_valid
    if (any(missing_rda)) {
      # Fit the same spatial smoother to display ranks if the saved values are absent.
      display_fit <- mgcv::gam(RDA_rank ~ s(x, y, bs = "tp", k = 120),
                              data = d[d$RDA_valid, ], method = "REML")
      d$RDA_display_rank[missing_rda] <- pmin(pmax(predict(display_fit,
        newdata = d[missing_rda, ]), 0), 1)
    }
    stopifnot(identical(d$RDA_display_rank[d$RDA_valid], d$RDA_rank[d$RDA_valid]),
      all(d$RDAforest_valid), all(is.finite(d$RDA_display_rank)))
    # Map full-domain percentile positions onto the shared-domain reference;
    # interpolated display values do not enter either rank calculation.
    on_shared_scale <- function(original, completed, reference) {
      ans <- reference
      ans[filled] <- approx(original[shared], reference[shared],
                            xout = completed[filled], rule = 2, ties = mean)$y
      stopifnot(all(is.finite(ans)), identical(ans[shared], reference[shared]))
      ans
    }
    linear <- on_shared_scale(d$RDA_rank, d$RDA_display_rank, d$linear_shared_rank)
    gf <- on_shared_scale(d$GF_rank, d$GF_display_rank, d$GF_shared_rank)
    rdf <- on_shared_scale(d$RDAforest_rank, d$RDAforest_rank, d$RDAforest_shared_rank)
    nonlinear <- (gf + rdf) / 2
    category <- rep("Intermediate", sum(filled))
    l <- linear[filled]; n <- nonlinear[filled]
    category[l >= .67 & n >= .67] <- "High–high"
    category[l <= .33 & n <= .33] <- "Low–low"
    category[l - n > .25] <- "Linear higher"
    category[n - l > .25] <- "Nonlinear higher"
    d$display_agreement <- d$agreement
    d$display_agreement[filled] <- category
    stopifnot(!anyNA(d$display_agreement),
      identical(d$display_agreement[shared], d$agreement[shared]))
    message("S6 plot only — GF interpolated: ", sum(!d$GF_valid),
      "; RDA interpolated for agreement: ", sum(missing_rda),
      "; RDAforest interpolated: 0; exploratory agreement: ", sum(filled),
      "; original agreement: ", sum(shared), "; novelty contour cells: ", sum(d$novel))
    p <- (offset_panel(d, "RDA", "Candidate-enriched RDA", "a", TRUE) |
          offset_panel(d, "GF", "Gradient Forest", "b", TRUE, completed = TRUE)) /
         (offset_panel(d, "RDAforest", "RDAforest", "c", TRUE) | agreement(d, completed = TRUE))
    export_plate(p, "Figure_S6_relative_genomic_offset_novelty_2071_2100_white_contour", 192)
  }
  for (relative in c(FALSE, TRUE)) {
    if (!((if (relative) "offset_rank" else "5") %in% tasks)) next
    for (d in list(mid, late)) {
      p <- (offset_panel(d, "RDA", "Candidate-enriched RDA", "a", relative) |
            offset_panel(d, "GF", "Gradient Forest", "b", relative)) /
           (offset_panel(d, "RDAforest", "RDAforest", "c", relative) | agreement(d))
      stem <- if (relative) {
        if (identical(d, mid)) "Figure_5_relative_genomic_offset_2041_2070" else
          "Figure_S6_relative_genomic_offset_novelty_2071_2100"
      } else {
        if (identical(d, mid)) "Figure_5_genomic_offset_2041_2070" else
          "Figure_S6_offset_novelty_2071_2100_number_pending"
      }
      export_plate(p, stem, 192)
    }
  }
}

# 5. Render supplementary structure and demographic figures ----

if ("S1" %in% tasks) {
  full <- read.csv(file.path(root, "04_PCA_14134_scores.csv"))
  rad <- read.csv(file.path(root, "04_PCA_one_per_RAD_scores.csv"))
  rad <- rad[match(full$InternalID, rad$InternalID), ]
  stopifnot(identical(full$InternalID, rad$InternalID), nrow(rad) == 245L)
  rad$ancestry_cluster <- full$ancestry_cluster
  for (pc in c("PC1", "PC2")) rad[[pc]] <- rad[[pc]] * sign(cor(rad[[pc]], full[[pc]]))
  panel <- function(x, title, tag) {
    ggplot(x, aes(PC1, PC2, colour = factor(ancestry_cluster))) + geom_point(size = 1.55) +
      scale_colour_manual(values = c(purple, "#B79AD2", "#F7C58D", "#C76100")) +
      labs(title = title, tag = tag, colour = "sNMF cluster") +
      theme_classic() + print_theme
  }
  p <- (panel(full, "Full panel: 14,134 SNPs", "a") |
        panel(rad, "One SNP per RAD locus", "b")) + plot_layout(guides = "collect")
  export_plate(p & theme(legend.position = "bottom"), "Figure_S1_PCA_RAD_representation", 96)
}

if ("demography" %in% tasks) {
  # Read TreeMix likelihoods and Stairway Plot summary tables.
  e <- new.env()
  expr <- parse("code/12_demography_treemix_figures.R")
  for (item in expr) {
    if (is.call(item) && as.character(item[[1]])[1] %in% c("dir.create", "write.csv", "ggsave")) next
    eval(item, e)
  }
  export_plate(e$p_llik + labs(subtitle = NULL) + print_theme, "Supplement_TreeMix_likelihood", 110)
  export_plate(e$p_stairway + print_theme, "Supplement_Stairway_Plot_2", 115)
}

if ("S7" %in% tasks) {
  s <- read.csv(file.path(root, "gradient_forest_all_14134/controlled_PC1_PC3_lon_lat/GF_controlled_predictor_importance_by_block.csv"))
  s$predictor <- reorder(toupper(s$predictor), s$importance, mean)
  p <- ggplot(s, aes(importance, predictor)) +
    geom_boxplot(fill = "#D8C9E5", colour = purple, width = .58) +
    stat_summary(fun = mean, geom = "point", colour = orange, orientation = "y") +
    labs(x = "Controlled GF predictor importance", y = NULL) + theme_classic() + print_theme
  export_plate(p, "Supplement_controlled_GF_importance_stability", 107)
}

# 6. Render the signed candidate-directionality figure ----

if ("directionality" %in% tasks) {
  Sys.setenv(AMEIVA_DEVELOPMENT_DIR = Sys.getenv("AMEIVA_DEVELOPMENT_DIR",
    "../Beyond-Generalism-development/candidate_geography_directionality"), AMEIVA_PROJECT_DIR = project)
  source_v5 <- file.path(Sys.getenv("AMEIVA_DEVELOPMENT_DIR"), "07_candidate_directionality_figure_v5.R")
  draw <- read_plot_code(source_v5, "pdf_file <-")
  paths <- file.path(out, paste0("Supplement_candidate_directionality", c(".pdf", ".png", ".tiff")))
  stopifnot(!any(file.exists(paths)))
  cairo_pdf(paths[1], width = 178 / 25.4, height = 234 / 25.4, bg = "white")
  draw(); dev.off()
  png(paths[2], width = 178, height = 234, units = "mm", res = 600, type = "cairo", bg = "white")
  draw(); dev.off()
  tiff(paths[3], width = 178, height = 234, units = "mm", res = 600, type = "cairo", bg = "white", compression = "lzw")
  draw(); dev.off()
}

# 7. Render structure and candidate diagnostics ----

if ("diagnostics" %in% tasks) {
  ent <- read.csv(file.path(root, "03_sNMF_cross_entropy_summary.csv"))
  p <- ggplot(ent, aes(K, mean)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "#D8C9E5") +
    geom_line(colour = purple) + geom_point(colour = purple) +
    geom_vline(xintercept = 4, linetype = 2, colour = orange) +
    scale_x_continuous(breaks = ent$K) + labs(x = "Number of ancestry components (K)",
      y = "Mean cross-entropy") + theme_classic() + print_theme
  export_plate(p, "Supplement_sNMF_cross_entropy", 100)

  model <- readRDS(file.path(root, "07_RDA_PCA_geo.rds"))
  perm <- readRDS(file.path(root, "07_RDA_PCA_geo_K4_nonquarter_geoonly_ecdf_v2_global_test.rds"))
  scan <- read.csv(file.path(root, "08_RDA_locus_scan_PCA_geo.csv"))
  stopifnot(length(attr(perm, "F.perm")) == 999, nrow(scan) == 14134)
  p_null <- ggplot(data.frame(F_perm = as.numeric(attr(perm, "F.perm"))), aes(F_perm)) +
    geom_histogram(bins = 32, fill = "#D8C9E5", colour = "white") +
    geom_vline(xintercept = perm$F[1], colour = orange) +
    labs(title = "Global permutation test", tag = "a", x = "Permuted F", y = "Permutations")
  bp <- as.data.frame(model$CCA$biplot[, 1:2])
  bp$predictor <- toupper(rownames(bp))
  p_load <- ggplot(bp, aes(RDA1, RDA2)) +
    geom_hline(yintercept = 0, colour = "grey88") + geom_vline(xintercept = 0, colour = "grey88") +
    geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), colour = purple,
      arrow = grid::arrow(length = unit(2, "mm"))) +
    geom_text(aes(label = predictor), colour = purple, nudge_y = .07, size = 3.2) +
    coord_equal(xlim = c(-1.13, 1.13), ylim = c(-1.13, 1.13)) +
    labs(title = "Climatic loadings", tag = "b")
  scan$detected <- ifelse(scan$candidate_Bonferroni01, "Bonferroni candidate", "Other SNP")
  scan$neglog <- -log10(pmax(scan$p_value, .Machine$double.xmin))
  p_detection <- ggplot(scan, aes(max_abs_loading_z, neglog, colour = detected)) +
    geom_point(size = .65, alpha = .55) +
    geom_hline(yintercept = -log10(unique(scan$bonferroni_threshold)), linetype = 2) +
    scale_colour_manual(values = c("Bonferroni candidate" = orange, "Other SNP" = purple), name = NULL) +
    labs(title = "Candidate detection", tag = "c", x = "Maximum absolute loading z-score", y = expression(-log[10](italic(P))))
  export_plate(((p_null | p_load) / p_detection) & theme_classic() & print_theme &
    theme(legend.position = "bottom"), "Supplement_pRDA_diagnostics", 176)

  folder <- file.path(root, "09_LFMM2_RDA_candidate_genotype_PCA")
  cand <- read.csv(file.path(folder, "03_RDA_LFMM2_shared_candidate_genotype_PCA_scores.csv"))
  k <- read.csv(file.path(folder, "03_RDA_LFMM2_shared_candidate_kmeans_diagnostic.csv"))
  size <- read.csv(file.path(folder, "00_candidate_set_PCA_kmeans_summary.csv"))
  stopifnot(identical(size$n_SNPs[size$candidate_set == "RDA_intersect_LFMM2"], 164L),
            nrow(cand) == 245L, all(cand$genotype_cluster %in% 1:2))
  p_scores <- ggplot(cand, aes(PC1, PC2, colour = factor(genotype_cluster))) +
    geom_point() + scale_colour_manual(values = c(purple, orange), name = "Genotype cluster") +
    labs(title = "Linear-consensus scores", tag = "a")
  p_k <- ggplot(k, aes(K, mean_silhouette)) + geom_line(colour = purple) +
    geom_point(colour = purple) + geom_point(data = subset(k, selected), colour = orange, size = 2) +
    scale_x_continuous(breaks = k$K) +
    labs(title = "Sensitivity to K", tag = "b", y = "Mean silhouette")
  export_plate((p_scores | p_k) & theme_classic() & print_theme & theme(legend.position = "bottom"),
    "Supplement_linear_consensus_K_diagnostics", 105)
}

# 8. Plot locus-specific allele-frequency profiles ----

if ("frequency" %in% tasks) {
  d <- read.csv(file.path(root, "11_164_RDA_LFMM2_locality_allele_frequency_environment.csv"))
  for (bio in c("bio2", "bio3", "bio13", "bio14")) {
    all <- subset(d, environmental_predictor == bio)
    loci <- sort(unique(all$panel_label))
    stopifnot(length(loci) > 0)
    groups <- split(loci, ceiling(seq_along(loci) / 15))
    for (i in seq_along(groups)) {
      q <- subset(all, panel_label %in% groups[[i]])
      p <- ggplot(q, aes(environmental_value, allele_frequency)) +
        geom_point(aes(size = n_individuals), colour = purple, alpha = .8) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = orange, linewidth = .55) +
        facet_wrap(~panel_label, ncol = 3) +
        scale_size_continuous(range = c(1.25, 3), name = "Individuals per locality") +
        scale_y_continuous(breaks = c(0,.5,1)) + coord_cartesian(ylim = c(0,1)) +
        labs(x = unique(q$environmental_variable), y = "Mean ALT frequency") +
        theme_minimal() + print_theme + theme(legend.position = "bottom",
          panel.grid.minor = element_blank(), strip.text = element_text(size = 9, colour = purple))
      export_plate(p, sprintf("Supplement_%s_ALT_frequency_part_%02d", toupper(bio), i),
        max(78, 36 * ceiling(length(groups[[i]]) / 3) + 35))
    }
  }
}

# 9. Plot offset magnitudes and RDAforest dimensional sensitivity ----

if ("offset_supplement" %in% tasks) {
  src <- "results/figures/RDA-2027-14134/Table_S_RDAforest_pcs8_pcs15_pixels.csv"
  d <- read.csv(src)
  panels <- c("Main: 8 axes" = "rank_pcs8", "Sensitivity: 15 axes" = "rank_pcs15",
              "Absolute rank difference" = "absolute_rank_difference")
  p <- do.call(rbind, lapply(names(panels), function(n)
    data.frame(d[c("x", "y", "period")], panel = n, value = d[[panels[[n]]]])))
  p$panel <- factor(p$panel, levels = names(panels))
  world <- map_data("world")
  map_bg <- geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE,
    fill = NA, colour = "grey70", linewidth = .18)
  xy <- coord_quickmap(xlim = c(-85,-30), ylim = c(-35,15))
  p_sens <- ggplot(p, aes(x, y, fill = value)) + geom_tile() + map_bg + xy +
    facet_grid(period ~ panel) +
    scale_fill_gradientn(colours = c(purple, "#F2EAF8", "#126E5F", "#FFF1E3", orange), limits = c(0, 1)) +
    labs(x = "Longitude", y = "Latitude", fill = "Rank / absolute difference") +
    theme_minimal() + print_theme + theme(panel.grid = element_blank(), legend.position = "bottom")
  export_plate(p_sens, "Supplement_RDAforest_8_15_axes", 137, width = 192)

  x <- read.csv("results/Figure_5_genomic_offset_final/Figure_5_relative_offset_source.csv")
  common <- subset(x, shared_valid)
  stopifnot(all(table(common$period) == c(3971, 2529)))
  cols <- c("2041-2070" = purple, "2071-2100" = orange)
  page <- lapply(c("RDA", "GF", "RDAforest"), function(method) {
    g <- data.frame(common[c("period", "x", "y")], offset = common[[method]])
    title <- if (method == "RDA") "Candidate-enriched RDA" else if (method == "GF") "Gradient Forest" else method
    ggplot(g, aes(offset, colour = period)) + stat_ecdf(linewidth = .7) +
      scale_colour_manual(values = cols, name = "Period") +
      labs(title = title, x = "Original genomic offset", y = "Cumulative fraction") +
      theme_classic() + print_theme
  })
  export_plate((wrap_plots(page, ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a")) & theme(legend.position = "bottom"),
    "Supplement_original_offset_distributions", 208)
}
