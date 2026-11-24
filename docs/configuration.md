# Paths and configuration

I have kept the validated analysis scripts unchanged. Most original scripts assign absolute paths inside their own source; the example below does **not** override those assignments. Re-running those scripts on another machine will require a separately reviewed, scientific-code-preserving path adaptation. [The path audit](path_audit.tsv) lists the affected files and lines.

The candidate-directionality v4 script reads four environment variables; its v5 figure-only script reads two. I can copy `config/paths.example.env` to a local untracked file, fill in my own paths, and source it before running either figure script. I must never commit the filled-in file.

```bash
cp config/paths.example.env .env
# Edit .env with local paths; this file is ignored by Git.
source .env
Rscript analysis/10_candidate_directionality/07_candidate_directionality_figure_v5.R
```

For v5, `AMEIVA_DEVELOPMENT_DIR` must contain the validated v4 summary tables, the external locality-frequency and locality-score tables, and `south_america_borders.gmt`. `AMEIVA_PROJECT_DIR` must contain `results/nonlinear_present_landscapes/GF_controlled_mean_cumulative_curves.csv`. The original v4 analysis additionally needs `AMEIVA_RESULTS_DIR` and `AMEIVA_ENV_FILE`. These variables are *input pointers*, not replacements for missing data permissions or provenance.

The full pipeline also requires the raw VCF, sample metadata, the 245 × 14,134 genotype matrices, CHELSA rasters, a land/range mask where applicable, saved GF/RDAforest objects, and genome assemblies/annotation. The current scripts reference both an external HD and another climate-raster directory. I do not infer GCM, SSP or download URLs from filenames; those fields remain open in the release checklist.
