# Gradient Forest methods and input-lineage record

The project contains two GF workflows with different response sets and climate predictors.

## A. Earlier candidate-locus GF and spatial offset

### Exact saved object and response set

The spatial maps and GF offset use:

```text
results/RDA-2027-14134/11_Gradient_Forest_nonquarter_geoonly_ecdf_v2.rds
```

The object contains 245 individuals, four environmental predictors and **118
candidate SNP responses**. Its response map is
`11_Gradient_Forest_SNP_map.csv` (118 rows). Of the 118 fitted responses, 108
have positive GF out-of-bag R² retained in the object `result` component; all
118 were nevertheless included as responses in the model fit.

The 118 SNPs were defined in `RDA-2027.R` as the intersection of:

1. FDR-significant outliers shared by the PCA + geography and sNMF + geography
   partial RDAs; and
2. candidates from the earlier LFMM analysis based on environmental PCs.

This earlier RDA ∩ LFMM consensus differs from the later 164-SNP RDA ∩ LFMM2 set and from the full-panel GF used for the Chapter 1 projections.

### Environmental predictors and preprocessing

The final saved object has exactly four predictors:

| Code | Environmental variable |
|---|---|
| bio2 | Mean diurnal range |
| bio3 | Isothermality |
| bio13 | Precipitation of wettest month |
| bio14 | Precipitation of driest month |

Raw CHELSA 1981--2010 values were matched to the 245 `InternalID`s and each
predictor was centred and scaled using its sample mean and sample standard
deviation. This was verified by reconstructing all four columns from
`Ameiva_Variables_rda.csv` and obtaining exact agreement with the saved
object's `X` matrix. Despite the file label `ecdf_v2`, **the fitted candidate
GF object contains z-standardised, not ECDF-transformed, predictors**.

Responses were the sNMF-imputed 0/1/2 genotype dosages. No additional GF-stage
residualisation for genomic PCs, ancestry coefficients, longitude or latitude
was applied. Population structure and geography entered the preceding partial RDA
candidate-selection stage, rather than this candidate-GF response matrix.

### GF parameters actually passed

The saved call in the RDS is:

```r
gradientForest(
  data = GF_data,
  predictor.vars = names(Environment),
  response.vars = GF_response_names,
  ntree = 1000,
  corr.threshold = 0.50,
  compact = TRUE,
  trace = TRUE,
  check.names = FALSE
)
```

`set.seed(2027)` was called immediately before the fit. `mtry`, `maxLevel`,
`nbin`, `transform`, and minimap-style options were not supplied in this call;
they therefore followed the installed `gradientForest` package defaults at the
time of the historical run. The original run did not save `sessionInfo()` or a
package-version manifest. A later project environment has
`gradientForest` 0.1-37, but the original run version remains unconfirmed.

### Spatial projection and domain handling

For the conservative original offset maps, CHELSA rasters were transformed by
the fitted cumulative-importance functions and cells outside the saved climate
domain were excluded. The separate full-coverage visualisation
`09b_GF_offset_full_endpoint_clamped.*` retains all terrestrial cells by
endpoint-clamping `cumimp` curves with `approx(..., rule = 2)`; its dashed teal
boundary marks novel climate and it is explicitly not a validated
extrapolation.

## B. Genome-wide raw and controlled GF (14,134-SNP analysis)

This separate workflow was run by `GF_14134_bruto_controlado_8cores.R` and is
stored under `gradient_forest_all_14134/`.

| Item | Recorded value |
|---|---:|
| Execution time | 2026-08-12 19:44:22 |
| Individuals | 245 |
| Imputed SNP responses | 14,134 |
| Represented RAD loci | 5,776 |
| Computational batches | 32 |
| Parallel workers used | 6 |
| Trees per response | 500 |
| `mtry` | 3 |
| `maxLevel` | 0 |
| `corr.threshold` | 0.50 |
| `compact` | TRUE |
| `nbin` | 101 |
| Seed | 2027 + block number; +10,000 for controlled blocks |

The ten predictors were `bio1`, `bio2`, `bio3`, `bio4`, `bio5`, `bio6`,
`bio7`, `bio12`, `bio13`, and `bio15`. They were transformed independently to
their empirical cumulative distribution functions using the 245 samples. The
raw model used imputed dosage responses. The controlled model first regressed
each dosage against PC1--PC3, standardised longitude and standardised latitude;
it then divided residuals by their SNP-specific residual standard deviation.
Both raw and controlled models were completed, and block objects were combined
after fitting.

This whole-panel GF serves a different purpose from the 118-locus candidate
GF: it is the genome-wide exploratory/nonlinear screening and relative-performance screening workflow. Its present and future projection sources are in `analysis/06_projection/`.

## Files inspected

- `RDA-2027.R`, GF and projection sections (earlier source)
- `11_Gradient_Forest_nonquarter_geoonly_ecdf_v2.rds`
- `11_Gradient_Forest_SNP_map.csv`
- `11_Gradient_Forest_importance.csv`
- `GF_14134_bruto_controlado_8cores.R`
- `gradient_forest_all_14134/00_GF_run_summary.csv`
