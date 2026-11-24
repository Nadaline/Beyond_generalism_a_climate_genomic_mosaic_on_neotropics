# Candidate climatic directionality and geography

The script `06_candidate_directionality_all_groups_v4.R` analyses eight signed LFMM2 response groups among the 164 SNPs shared by partial RDA and LFMM2. The companion script `07_candidate_directionality_figure_v5.R` draws the signed associations, observed frequencies, maps and GF turnover curves from the saved summaries.

## Candidate directionality figure

The earlier v5 export is 176 × 234 mm; its caption is in [the figure notes](../../docs/candidate_directionality_v5/). The Chapter 1 figure is [Figure 3](../../results/figures/Figure_3_candidate_directionality.pdf). The exact correspondence of the two selected map groups is discussed in [limitations](../../docs/known_limitations.md).

| Panel | Saved input | Interpretation |
| --- | --- | --- |
| a | `30_v4_group_counts.csv`, `31_v4_group_mean_signed_effects.csv` | Mean signed LFMM2 ALT-dosage effects for all eight groups and four predictors; row counts are SNPs/distinct RAD loci. |
| b | External `32_v4_all_group_locality_ALT_frequency.csv`, `33_v4_observed_frequency_by_climate_decile.csv` | Locality observations and ten individual-count-weighted descriptive climate-quantile means, not fitted LFMM2 curves. |
| c | External `18_LFMM2_group_locality_scores_v2.csv`, `38_v4_map_display_scale.csv`, `39_v4_main_map_group_selection.csv`, `south_america_borders.gmt` | Two selected observed-site maps, within-group standardised and not spatially interpolated. |
| d | External `results/nonlinear_present_landscapes/GF_controlled_mean_cumulative_curves.csv` | Curves from the **controlled 14,134-SNP GF**, not the 164 linear candidates or the historical 118-SNP GF. |

The repository contains the small tables 30, 31, 33, 38 and 39, the border file and the v5 script, but not the exact-location tables 18/32 or the GF curve file. Full regeneration requires access to those inputs (see [data availability](../../docs/data_availability.md)).

## Running the figure-only script

Set `AMEIVA_DEVELOPMENT_DIR` to the directory holding the v4 tables and the border file, and `AMEIVA_PROJECT_DIR` to the project root containing the saved GF curve file. `config/paths.example.env` shows the variables without exposing a personal path. After sourcing a filled-in local `.env`, run:

```bash
Rscript analysis/10_candidate_directionality/07_candidate_directionality_figure_v5.R
```

The script checks expected input dimensions and the two map groups before drawing. It needs `ggplot2`, `scales`, a Cairo-capable R graphics device and base `grid`. It writes v5 PDF/PNG into `AMEIVA_DEVELOPMENT_DIR`; it does not write to `figures/` automatically.

## Original v4 analysis

The v4 script additionally reads `AMEIVA_RESULTS_DIR` (the 14,134-SNP results) and `AMEIVA_ENV_FILE` (the individual CHELSA 1981–2010 table). It verifies 245 aligned individuals, the partial RDA–LFMM2 intersection and original VCF ALT coding. Group frequencies are means of imputed dosage/2 across SNPs and individuals at each coordinate-defined locality. It uses fractional-response quasi-binomial GLMs for a conditional diagnostic and existing checkerboard folds for sensitivity, not independent confirmation of selection. The original candidate set and genome-wide PCs preceded those folds. A one-SNP-per-RAD sensitivity is saved separately.

BIO2 is displayed in °C, BIO3 as a ratio, and BIO13/BIO14 in mm. ALT is reference-relative, not a beneficial allele. A marginal observed frequency trend and a conditional LFMM2 coefficient answer different questions. The GF cumulative curves indicate turnover importance, not the direction of allele-frequency change or a fitness effect.
