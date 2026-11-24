# Software and dependency record

The table distinguishes declared dependencies from versions recorded for individual analyses. A complete environment record is not available.

| Component | Script-declared dependencies | Verified version information |
| --- | --- | --- |
| Core filtering, structure, partial RDA, LFMM2, RDAforest and annotation | R; `vcfR`, `dartR.base`, `dartR.popgen`, `dartR.spatial`, `adegenet`, `vegan`, `LEA`, `Boruta`, `usdm`, `qvalue`, `robust`, `ggplot2`, `cluster`, `terra`, `gradientForest`, `RDAforest`, `GenomicRanges`, `IRanges`, `maps`, `sf`; external `bcftools` | The LFMM2 signed-effect recovery file records R 4.4.0 and LEA 3.22.0 for that recovery. Versions for other steps remain unconfirmed. |
| Full-panel GF | R; `gradientForest`, `parallel`, `LEA`, `ggplot2`; later plotting sections also load `ggVennDiagram`, `VennDiagram`, `grid` | The GF run-summary records parameters and completion, not an R/package version lockfile. |
| RDAforest | R; `RDAforest`, `vegan` and dependencies of the core script | A local `rdaforest_session_info.txt` records R 4.5.3 and RDAforest 2.6.10, although its correspondence to the exact main fit is unconfirmed. |
| Climate projection and maps | R; `terra`, `parallel`, `ggplot2`, `patchwork`, `maps`, `gradientForest`; Figure 4's external finishing script also loads `mgcv` | Original versions and climate-source metadata still need a run-level manifest. |
| Liftoff | Liftoff, Python environment, `curl`, `gzip`; annotation in R uses `GenomicRanges` and `IRanges` | The output GFF records Liftoff 1.6.3. The transfer command is preserved in `provenance/LIFTOFF_METHODS_RECORD.md`. |
| Figure 2 v5 | R; `ggplot2`, `scales`, base `grid`, Cairo-capable graphics device | The v4 README records R 4.4.0, ggplot2 4.0.3 and scales 1.4.0 for its figure environment; the v5 export was checked for embedded PDF fonts but has no saved session manifest. |

The core script also invokes or prepares external tools for TreeMix and Stairway Plot 2; their executable versions and full invocation environments still need to be recorded before public release. 
