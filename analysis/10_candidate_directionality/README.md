# Candidate climatic directionality and geography

I use `06_candidate_directionality_all_groups_v4.R` for the scientific analysis of eight signed LFMM2 response groups among the 164 SNPs shared by the final partial-RDA and LFMM2 candidate sets. I retain its v4 main and supplementary figures and summary tables as validated development history. I use `07_candidate_directionality_figure_v5.R` **only to redraw** the current intended manuscript Figure 2 from those saved tables; it does not rerun LFMM2, pRDA, GF, spatial GLMs or candidate selection.

## Current Figure 2

The current private figure is `figures/candidate_directionality_v5/Figure_candidate_directionality_v5.pdf` (also PNG), with its caption in `docs/candidate_directionality_v5/`. I exported it at 176 × 234 mm. The v4 PDF/PNG and caption remain under their v4 paths and are not silently replaced.

| Panel | Saved input | Interpretation |
| --- | --- | --- |
| a | `30_v4_group_counts.csv`, `31_v4_group_mean_signed_effects.csv` | Mean signed LFMM2 ALT-dosage effects for all eight groups and four predictors; row counts are SNPs/distinct RAD loci. |
| b | External `32_v4_all_group_locality_ALT_frequency.csv`, `33_v4_observed_frequency_by_climate_decile.csv` | Locality observations and ten individual-count-weighted descriptive climate-quantile means, not fitted LFMM2 curves. |
| c | External `18_LFMM2_group_locality_scores_v2.csv`, `38_v4_map_display_scale.csv`, `39_v4_main_map_group_selection.csv`, `south_america_borders.gmt` | Two selected observed-site maps, within-group standardized and not spatially interpolated. |
| d | External `results/nonlinear_present_landscapes/GF_controlled_mean_cumulative_curves.csv` | Curves from the **controlled 14,134-SNP GF**, not the 164 linear candidates or the historical 118-SNP GF. |

The repository contains the small tables 30, 31, 33, 38 and 39, the border file and the v5 script, but not the exact-location tables 18/32 or the GF curve file. A third party therefore cannot fully regenerate v5 from this checkout alone. I will decide an access route for these inputs before public release.

## Running the figure-only script

I set `AMEIVA_DEVELOPMENT_DIR` to the directory holding all validated v4 tables and the border file, and `AMEIVA_PROJECT_DIR` to the project root containing the saved GF curve file. `config/paths.example.env` shows the variables without exposing a personal path. After sourcing a filled-in local `.env`, I can run:

```bash
Rscript analysis/10_candidate_directionality/07_candidate_directionality_figure_v5.R
```

The script checks expected input dimensions and the two map groups before drawing. It needs `ggplot2`, `scales`, a Cairo-capable R graphics device and base `grid`. It writes v5 PDF/PNG into `AMEIVA_DEVELOPMENT_DIR`; it does not write to `figures/` automatically.

## Original v4 analysis

The v4 script additionally reads `AMEIVA_RESULTS_DIR` (the validated 14,134-SNP results) and `AMEIVA_ENV_FILE` (the individual CHELSA 1981–2010 table). It verifies 245 aligned individuals, the current pRDA–LFMM2 intersection and original VCF ALT coding. Group frequencies are means of imputed dosage/2 across SNPs and individuals at each coordinate-defined locality. It uses fractional-response quasi-binomial GLMs for a conditional diagnostic and existing checkerboard folds for sensitivity, not independent confirmation of selection. The original candidate set and genome-wide PCs preceded those folds. A one-SNP-per-RAD sensitivity is saved separately.

BIO2 is displayed in °C, BIO3 as a ratio, and BIO13/BIO14 in mm. ALT is reference-relative, not a beneficial allele. A marginal observed frequency trend and a conditional LFMM2 coefficient answer different questions. The GF cumulative curves indicate turnover importance, not the direction of allele-frequency change or a fitness effect.
