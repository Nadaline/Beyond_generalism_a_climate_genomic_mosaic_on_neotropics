# Interpreting the saved figures and models

The six figure previews come from the R7 thesis. Their PDF companions are viewing copies made from those images, rather than the original plotting exports. The scientific scripts are available in `analysis/` and `scripts/`, but saved model objects, original RADseq inputs and some environmental rasters are not all part of the repository.

The projected genomic offsets represent method-specific positions, not fitness predictions. The R7 map uses relative ranks on a 0–1 scale; direct comparison of raw distances between methods is inappropriate. Earlier figure checks also identified differences between a reported spatial correlation and a saved comparison grid, and display-only filling in the late-century supplementary map. These issues warrant reconciliation with the original fitted objects, climate masks and comparison grid before reuse of those estimates. Figure rendering and model fitting are separate operations.

## Additional reproducibility checks

- The public renderer refers to `scripts/make_Figure3_nonlinear_models_v4.R` and `scripts/finalize_Figure4_present_nonlinear_landscapes.R`; these sources are available under `analysis/09_figures/`. Its `code/12_demography_treemix_figures.R` dependency is not deposited. The original executable paths are retained.
- The signed-group map-selection table and v5 script identify BIO14+ and BIO13−. The previous README named BIO3− instead; the exact Chapter 1 figure-to-input correspondence remains to be checked.
- The late-century GF finishing script fills unsupported cells for display, but its generated caption describes only supported cells. Interpret the saved map with the domain mask and the [projection record](figure_5_input_output.md).
- SVDquartets and the complete TreeMix/Stairway Plot execution workflows are not available as standalone scripts in this checkout. The core script contains input preparation for TreeMix and Stairway Plot 2.

For reproducibility at each level, see the [reproduction guide](reproducibility.md) and [technical input record](technical_reproduction_notes.md).
