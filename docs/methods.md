# From RAD markers to climatic landscapes

## Which response is being modelled?

The principal marker matrix contains 245 individuals, 14,134 SNPs and 5,776 RAD loci. The K = 4 genetic clusters describe broad population structure, not the ten geographic groups used in TreeMix. The partial RDAs and LFMM2 use **allelic dosage as the response**, with climate predictors BIO2, BIO3, BIO13 and BIO14. The two geographically conditioned partial RDAs use different structure covariates (genomic PCs or sNMF proportions). Their loading scans share 889 SNPs; the stricter subset has 626. The LFMM2 result is the union of per-predictor FDR discoveries (314 SNPs); the intersection with the 889-SNP partial RDA list is 164. See the [conditioned fits](../results/summary_tables/Table_1_R7_partial_RDA.csv) and [candidate counts](../results/summary_tables/R7_candidate_counts_verified.csv).

## How do the nonlinear models differ?

The [GF script](../analysis/04_gradient_forest/GF_14134_bruto_controlado_8cores.R) distributes 14,134 markers across 32 batches and fits 500 trees per marker at `mtry = 3`. The controlled response removes PC1–PC3 and geographic coordinates before scaling; the ten environmental predictors use empirical percentiles. Variable importance and the cumulative turnover curves describe the fitted relationships. The 707 SNPs in the upper 5% of out-of-bag R² form a relative-performance subset, not a significance-based set of independent discoveries. The 20 GF-only loci were a separate exploratory shortlist; four retain estimable response variation in all four checkerboard partitions. That is not evidence of adaptive function.

RDAforest represents genome-wide covariance with genomic ordination axes; its main fit uses eight axes and retained BIO2, BIO3 and BIO13. The 15-axis fit is a sensitivity comparison. Present-day PC1 and future offset are different summaries. RDAforest's offset uses its own fitted transformation and present q90 scale, not the GF cumulative-turnover distance or the eigenvalue-weighted RDA distance.

## What do annotations add?

[Liftoff](../analysis/07_annotation/prepare_liftoff_ameiva.sh) transfers features between assemblies; the core R source compares intervals with candidate SNPs. A reported overlap or nearest annotated gene gives genomic context, not functional confirmation. For a different genome verify that SNP and feature coordinates use the same assembly and sequence names before calculating overlaps.
