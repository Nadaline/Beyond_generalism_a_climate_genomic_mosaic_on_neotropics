# Working from figures, processed inputs or reads

## A. Inspect or redraw figures from saved values

The [figure index](figures.md) identifies the files read by each panel. To redraw Figure 1 in the complete project tree, supply the authorised PCA model and scores, sNMF proportions and environment table, set `AMEIVA_RESULTS_DIR` and `AMEIVA_ENV_FILE`, create an empty output directory and set `AMEIVA_PUBLICATION_DIR` to it. Run `Rscript scripts/render_ameiva_publication_figures.R 1` from that tree to write the population-structure PDF, PNG and TIFF. For Figure 2 use task `2` with vegan, GDAL, the saved conditioned ordination, candidate PCA and raster. The other tasks and companion plot files are listed in the script. Some companion paths require adaptation (see [limitations](known_limitations.md)). Only reference exports and some summary tables are included here; these redraws cannot be run from a fresh checkout alone.

## B. Refit from processed genotypes

Obtain the ordered dosage matrix or filtered Stacks VCF, RAD identifiers, locality and climate records, and structure covariates (see [data access](data_availability.md)). Start with the [core analysis](../analysis/01_core/RDA-2027-bonferroni.R), verify 245 unique ordered samples, 14,134 marker IDs and 5,776 RAD loci, then use its structure, partial-RDA and LFMM2 sections before calculating intersections. Run the [full-panel GF](../analysis/04_gradient_forest/GF_14134_bruto_controlado_8cores.R) separately; retain the 32 fitted batches as well as their OOB and importance tables. For spatial projections, also retain the fit, ordered climate layers, comparison grid and masks. The exact R7 projection inputs are not all available here.

Software is method-specific: R with vegan, LEA and the declared core packages; gradientForest and RDAforest for the nonlinear models; terra, ggplot2, patchwork and GDAL for maps. The LFMM2 effects recovery records R 4.4.0 and LEA 3.22.0, and the annotation GFF records Liftoff 1.6.3; these observations do not pin the environments of every run. See the [dependency record](../docs/dependencies.md). Save `sessionInfo()` alongside new fitted objects instead of assuming a version for earlier runs.

## C. Reconstruct from reads

The core script starts at a Stacks VCF. Rebuilding that input requires the read accessions and barcodes, the reference, and the demultiplexing, alignment and variant-filtering settings. This repository does not contain a verified reads-to-VCF recipe or all authorised inputs. The earlier scripts describe individual stages rather than a complete read-to-figure workflow.
