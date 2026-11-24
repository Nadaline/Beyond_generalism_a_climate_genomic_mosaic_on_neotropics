# Interpreting the saved figures and models

The six figure previews come from the R7 thesis. Their PDF companions are viewing copies made from those images, rather than the original plotting exports. The scientific scripts are available in `analysis/` and `scripts/`, but saved model objects, original RADseq inputs and some environmental rasters are not all part of the repository.

The projected genomic offsets represent method-specific positions, not fitness predictions. The R7 map uses relative ranks on a 0–1 scale; direct comparison of raw distances between methods is inappropriate. Earlier figure checks also identified differences between a reported spatial correlation and a saved comparison grid, and display-only filling in the late-century supplementary map. These issues warrant reconciliation with the original fitted objects, climate masks and comparison grid before reuse of those estimates. The original models should not be refitted merely to change how the maps look.

For reproducibility at each level, see the [reproduction guide](reproducibility.md) and [technical input record](technical_reproduction_notes.md).
