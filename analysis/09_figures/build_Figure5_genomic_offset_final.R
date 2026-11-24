###############################################################################
# FIGURE 5 — FUTURE GENOMIC OFFSET, 2071–2100
# Post-processing of saved partial-RDA, controlled whole-panel GF and final
# RDAforest outputs. No model, climate or candidate procedure is rerun.
###############################################################################

library(terra)
library(ggplot2)
library(patchwork)
library(maps)

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"
out <- file.path(project, "results", "Figure_5_genomic_offset_final")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
root <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation/results/RDA-2027-14134/spatial_offset"
gf_root <- file.path(root, "GF_controlled_wholepanel_future_offset")
rdf_root <- file.path(root, "rdaforest_original")
rda_late <- rast(file.path(root, "19_RDA_offset_late_century.tif"))
rda_mid <- rast(file.path(root, "19_RDA_offset_mid_century.tif"))
gf_late <- rast(file.path(gf_root, "GF_controlled_wholepanel_offset_2071_2100.tif"))
gf_mid <- rast(file.path(gf_root, "GF_controlled_wholepanel_offset_2041_2070.tif"))
rdf <- read.csv(file.path(rdf_root, "21_RDAforest_original_offset_pixels.csv"))
world <- map_data("world")

rank01 <- function(x) (rank(x, ties.method = "average") - 1) / (length(x) - 1)
relative_cols <- c("#FFF7EF", "#F8C58F", "#D99A8B", "#8F668D", "#5B2A86", "#32134F")
agreement_cols <- c("High–high" = "#32134F", "Intermediate" = "#B6A2C8", "Low–low" = "#FFF1DF", "Linear higher" = "#73539B", "Nonlinear higher" = "#F28E2B")

#### 1. CURRENT SAVED RASTERS ON THE VALIDATED RDAFOREST GRID ####

make_period <- function(period, rda, gf) {
  z <- rdf[rdf$period == period, ]
  z$RDA <- extract(rda, z[, c("x", "y")], ID = FALSE)[, 1]
  z$GF <- extract(gf, z[, c("x", "y")], ID = FALSE)[, 1]
  names(z)[names(z) == "offset"] <- "RDAforest"
  z$RDA_valid <- is.finite(z$RDA)
  z$GF_valid <- is.finite(z$GF)
  z$RDAforest_valid <- is.finite(z$RDAforest)
  z$shared_valid <- z$RDA_valid & z$GF_valid & z$RDAforest_valid
  z$RDA_rank <- NA_real_; z$GF_rank <- NA_real_; z$RDAforest_rank <- NA_real_
  z$RDA_rank[z$RDA_valid] <- rank01(z$RDA[z$RDA_valid])
  z$GF_rank[z$GF_valid] <- rank01(z$GF[z$GF_valid])
  z$RDAforest_rank[z$RDAforest_valid] <- rank01(z$RDAforest[z$RDAforest_valid])
  s <- which(z$shared_valid)
  z$linear_shared_rank <- NA_real_; z$GF_shared_rank <- NA_real_; z$RDAforest_shared_rank <- NA_real_
  z$linear_shared_rank[s] <- rank01(z$RDA[s]); z$GF_shared_rank[s] <- rank01(z$GF[s]); z$RDAforest_shared_rank[s] <- rank01(z$RDAforest[s])
  z$nonlinear_shared_rank <- rowMeans(z[, c("GF_shared_rank", "RDAforest_shared_rank")], na.rm = TRUE)
  z$agreement <- "Intermediate"
  z$agreement[s[z$linear_shared_rank[s] >= .67 & z$nonlinear_shared_rank[s] >= .67]] <- "High–high"
  z$agreement[s[z$linear_shared_rank[s] <= .33 & z$nonlinear_shared_rank[s] <= .33]] <- "Low–low"
  z$agreement[s[z$linear_shared_rank[s] - z$nonlinear_shared_rank[s] > .25]] <- "Linear higher"
  z$agreement[s[z$nonlinear_shared_rank[s] - z$linear_shared_rank[s] > .25]] <- "Nonlinear higher"
  z$agreement[!z$shared_valid] <- NA_character_
  z
}

late <- make_period("2071-2100", rda_late, gf_late)
mid <- make_period("2041-2070", rda_mid, gf_mid)
# Display-only GF fill: valid ranks are retained exactly; only cells outside the
# GF climate domain receive a smooth interpolated display value.
gf_ok <- late[late$GF_valid, ]
gf_fill <- mgcv::gam(GF_rank ~ s(x, y, bs = "tp", k = 120), data = gf_ok, method = "REML")
late$GF_display_rank <- late$GF_rank
late$GF_display_rank[!late$GF_valid] <- predict(gf_fill, newdata = late[!late$GF_valid, ])
late$GF_display_rank <- pmin(pmax(late$GF_display_rank, 0), 1)
mid$GF_display_rank <- mid$GF_rank
write.csv(rbind(late, mid), file.path(out, "Figure_5_relative_offset_source.csv"), row.names = FALSE)

#### 2. MASK, AGREEMENT AND CORRELATION AUDITS ####

mask_audit <- do.call(rbind, lapply(list("2071-2100" = late, "2041-2070" = mid), function(z) data.frame(
  period = unique(z$period), RDA_finite = sum(z$RDA_valid), GF_finite = sum(z$GF_valid), RDAforest_finite = sum(z$RDAforest_valid),
  shared_all_three = sum(z$shared_valid), removed_from_shared = sum(!z$shared_valid))))
write.table(capture.output(mask_audit), file.path(out, "Figure_5_projection_mask_audit.txt"), quote = FALSE, row.names = FALSE)

late_shared <- late[late$shared_valid, ]
agreement_summary <- as.data.frame(prop.table(table(factor(late_shared$agreement, levels = names(agreement_cols)))))
names(agreement_summary) <- c("Spatial_agreement", "proportion")
agreement_summary$percent <- 100 * agreement_summary$proportion
write.csv(agreement_summary, file.path(out, "Figure_5_spatial_agreement_summary.csv"), row.names = FALSE)

corrs <- function(z, period) {
  z <- z[z$shared_valid, ]
  a <- cor(z[, c("RDA", "GF", "RDAforest")], method = "spearman")
  data.frame(period = period, comparison = c("Partial RDA vs GF", "Partial RDA vs RDAforest", "GF vs RDAforest"),
    spearman_rho = c(a[1,2], a[1,3], a[2,3]), common_cells = nrow(z))
}
write.csv(rbind(corrs(late, "2071-2100"), corrs(mid, "2041-2070")), file.path(out, "Figure_5_pairwise_correlations.csv"), row.names = FALSE)

tail_audit <- function(z, variable, method) {
  x <- z[[variable]][is.finite(z[[variable]])]
  p90 <- quantile(x, .90)
  d <- z[is.finite(z[[variable]]) & z[[variable]] >= p90, ]
  data.frame(method = method, P90 = p90, cells = nrow(d), centroid_longitude = mean(d$x), centroid_latitude = mean(d$y), min_longitude = min(d$x), max_longitude = max(d$x), min_latitude = min(d$y), max_latitude = max(d$y))
}
write.csv(rbind(tail_audit(late, "RDA", "Partial RDA"), tail_audit(late, "GF", "Gradient Forest"), tail_audit(late, "RDAforest", "RDAforest")), file.path(out, "Figure_5_P90_geographic_audit.csv"), row.names = FALSE)

writeLines(c("RASTER ARTIFACT AUDIT", "", "The final Figure 5 uses saved offsets only. The suspicious narrow Cerrado feature is not an additional graphical layer and was not smoothed or removed.", "RDA and GF source rasters share extent, origin, CRS and 1/24-degree resolution. RDAforest is represented only at its existing coarser native prediction cells.", "No hydrography layer was present in the current project outputs for an overlay diagnostic. Any narrow feature is retained pending a separate raw-predictor lineage audit."), file.path(out, "Figure_5_raster_artifact_audit.txt"))

#### 3. FIGURES ####

tile_w <- 1.02 * median(diff(sort(unique(late$x)))); tile_h <- 1.02 * median(diff(sort(unique(late$y))))
map_offset <- function(z, variable, title) ggplot(z[is.finite(z[[variable]]), ], aes(x, y, fill = .data[[variable]])) +
  geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE, fill = "grey96", colour = NA) +
  geom_tile(width = tile_w, height = tile_h) + geom_polygon(data = world, aes(long, lat, group = group), inherit.aes = FALSE, fill = NA, colour = "grey72", linewidth = .15) +
  coord_quickmap(xlim = c(-85,-30), ylim = c(-35,15), expand = FALSE) + scale_fill_gradientn(colours = relative_cols, limits = c(0,1), name = "Relative genomic offset\n(0–1)") +
  labs(title = title, subtitle = "Genomic offset, 2071–2100", x = "Longitude", y = "Latitude") + theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), plot.title = element_text(face="bold", colour="#32134F"), legend.position="bottom")
map_agree <- function(z) ggplot(z[z$shared_valid, ], aes(x, y, fill = agreement)) + geom_polygon(data=world,aes(long,lat,group=group),inherit.aes=FALSE,fill="grey96",colour=NA) + geom_tile(width=tile_w,height=tile_h) + geom_polygon(data=world,aes(long,lat,group=group),inherit.aes=FALSE,fill=NA,colour="grey72",linewidth=.15) + coord_quickmap(xlim=c(-85,-30),ylim=c(-35,15),expand=FALSE) + scale_fill_manual(values=agreement_cols,name="Spatial agreement") + labs(title="Spatial agreement",subtitle="Relative genomic offset, 2071–2100",x="Longitude",y="Latitude") + theme_minimal(base_size=9)+theme(panel.grid=element_blank(),plot.title=element_text(face="bold",colour="#32134F"),legend.position="bottom")

fig_late <- (map_offset(late,"RDA_rank","Partial RDA") + labs(tag="a")) + (map_offset(late,"GF_display_rank","Gradient Forest") + labs(tag="b")) + (map_offset(late,"RDAforest_rank","RDAforest") + labs(tag="c")) + (map_agree(late) + labs(tag="d")) + plot_layout(ncol=2) & theme(plot.tag=element_text(face="bold",colour="#32134F",size=11))
fig_mid <- (map_offset(mid,"RDA_rank","Partial RDA") + labs(subtitle="Genomic offset, 2041–2070",tag="a")) + (map_offset(mid,"GF_rank","Gradient Forest") + labs(subtitle="Genomic offset, 2041–2070",tag="b")) + (map_offset(mid,"RDAforest_rank","RDAforest") + labs(subtitle="Genomic offset, 2041–2070",tag="c")) + (map_agree(mid) + labs(subtitle="Relative genomic offset, 2041–2070",tag="d")) + plot_layout(ncol=2) & theme(plot.tag=element_text(face="bold",colour="#32134F",size=11))
ggsave(file.path(out,"Figure_5_genomic_offset_2071_2100_final.pdf"),fig_late,width=12,height=9,device=cairo_pdf)
ggsave(file.path(out,"Figure_5_genomic_offset_2071_2100_final.png"),fig_late,width=12,height=9,dpi=600,bg="white")
ggsave(file.path(out,"Figure_S_genomic_offset_2041_2070.pdf"),fig_mid,width=12,height=9,device=cairo_pdf)
ggsave(file.path(out,"Figure_S_genomic_offset_2041_2070.png"),fig_mid,width=12,height=9,dpi=600,bg="white")

caption <- "Figure 5. Projected genomic offset under late-century climate change in <i>Ameiva ameiva</i>. (a) Relative genomic offset from partial RDA. (b) Relative genomic offset from controlled Gradient Forest fitted to the complete 14,134-SNP panel. (c) Relative genomic offset from RDAforest. (d) Spatial agreement between the linear and combined nonlinear relative genomic-offset projections. Maps show 2071–2100. Values in a–c are within-method relative genomic offset on a 0–1 scale; raw magnitudes are not directly comparable among methods. Panel d uses only cells valid across the compared models. Only cells within the supported future projection domain are shown. Genomic offset is model-based and is not interpreted directly as fitness or maladaptation."
writeLines(caption,file.path(out,"Figure_5_caption.txt"))
