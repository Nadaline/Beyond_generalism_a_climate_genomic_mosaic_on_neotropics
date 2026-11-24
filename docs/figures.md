# Figures in Chapter 1

The R7 numbering is used here even when an earlier filename has an earlier number. The PNG and PDF exports for Figures 1–6 are the chapter-matched versions; similarly named earlier scripts and exports are not interchangeable.

| R7 figure | PDF | Script and the inputs that determine the panels |
|---|---|---|
| 1. Continental population structure | [PDF](../results/figures/Figure_1_continental_structure.pdf) | [Task 1](../scripts/render_ameiva_publication_figures.R) uses the [plotting source](../scripts/restyle_main_figures_1_2.R), saved PCA scores and model, sNMF proportions and sample geography. |
| 2. Structure-controlled linear variation | [PDF](../results/figures/Figure_2_linear_climate_associations.pdf) | [Task 2](../scripts/render_ameiva_publication_figures.R) combines the saved `07_RDA_PCA_geo.rds` (a), PCA and cluster scores for the 164 shared SNPs (b–c), and two bands of `18_RDA_present_axes.tif` (d–e). The raster is averaged in 4 × 4 blocks for display. |
| 3. Candidate associations and geography | [PDF](../results/figures/Figure_3_candidate_directionality.pdf) | [v5 plot source](../analysis/10_candidate_directionality/07_candidate_directionality_figure_v5.R). The inputs distinguish signed LFMM2 effects, ALT frequencies, local genotype profiles and the separately fitted GF turnover. |
| 4. Nonlinear genomic mosaic | [PDF](../results/figures/Figure_4_nonlinear_mosaic.pdf) | [Task 3](../scripts/render_ameiva_publication_figures.R) reuses the [earlier panel script](../analysis/09_figures/make_Figure3_nonlinear_models_v4.R): GF batch importance, 707 turnover summaries, four GF-only examples and RDAforest jackknife values. |
| 5. Present-day nonlinear landscapes | [PDF](../results/figures/Figure_5_present_nonlinear_landscapes.pdf) | [Task 4](../scripts/render_ameiva_publication_figures.R) uses the [present-day map and inset source](../analysis/09_figures/finalize_Figure4_present_nonlinear_landscapes.R), native PC1 rasters and the saved climatic correlations. RDAforest bilinear resampling is for display only. |
| 6. Mid-century genomic offset | [PDF](../results/figures/Figure_6_relative_genomic_offset_2041_2070.pdf) | The [offset plot source](../analysis/09_figures/build_Figure5_genomic_offset_final.R) traces the three saved offset surfaces into a comparison grid. The later [`offset_rank` task](../scripts/render_ameiva_publication_figures.R) reads the saved comparison grid for the manuscript plate. |

## Tracing the RDA surfaces

The [core script, sections 18–19](../analysis/01_core/RDA-2027-bonferroni.R), documents this chain:

```text
PC/geography and sNMF/geography pRDAs → 889 shared candidates
strict loading criterion → 626 SNPs → climate-only RDA
  ├─ 18_RDA_projection_model.rds → 18_RDA_present_axes.tif → Figure 2d–e
  └─ weighted present/future distances → 19_RDA_offset_mid_century.tif
     → saved model-comparison grid → mid-century Figure 6
```

The code explicitly subsets `candidates_RDA_material_strict` in that fit; neither map is a direct projection of the all-SNP conditioned ordination. The R7 describes the 889 candidates as the landscape-RDA panel. The script path points to 626 for **both** projected products, but the final renderers read named rasters and a comparison table, not the fit itself. Confirming the final inputs requires the exact `18_RDA_projection_model.rds`, its response SNP IDs, and the matching present raster, future offset raster and saved comparison grid to distinguish the documented command from any later replacement. Compare the fit's response IDs with `08_RDA_candidates_landscape_consensus.csv` and `08_RDA_candidates_landscape_strict_3SD.csv`; link the raster and grid checksums to those used in the plate.

## Projected offset

The Figure 6 colour bars are method-specific ranks on 0–1; the R7 caption instead describes original distances. R7 reports mid-century RDA–GF Spearman ρ = 0.675; the inspected saved grid yielded 0.7270 on 3,971 shared valid cells. We do not substitute that value for the R7 statistic without identifying the exact grid and masks used. Figure S5 (2071–2100) uses display-only filling in parts of the GF and agreement maps; the contour does not make those cells model-supported. Raw magnitudes are a separate within-method comparison (Figure S9). Unsupported and missing cells have different meanings in these plates. The [earlier offset notes](../docs/figure_5_input_output.md) describe differences among the saved surfaces.
