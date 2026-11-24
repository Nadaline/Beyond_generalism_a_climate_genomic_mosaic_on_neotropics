# Data availability

The code, figures and selected aggregate tables are available here. Other computational inputs remain outside this repository because they require a separate deposit, access route or verified acquisition record. An external file named in a script is not automatically an openly deposited dataset.

| Input or artifact | Current role | Access or provenance information still needed |
| --- | --- | --- |
| Stacks RADseq VCF `Ameiva_genome/populations.snps.vcf`, raw reads and sample metadata | Original genotype calls and the validated 245-individual subset | An accession or approved archive and a sample-identifier review are needed. |
| Imputed 245 × 14,134 genotype matrix, `.lfmm`, saved PCA/sNMF objects and candidate tables | Inputs for partial RDA, LFMM2, GF, RDAforest and figure regeneration | A versioned deposit with checksums and column order is needed for these large or individual-level inputs. |
| CHELSA present and future raster layers | BIO predictors for model fitting and projections | Required metadata include the exact CHELSA product, original resolution, GCM, SSP, acquisition links, licenses and local aggregation. The present scripts alone do not establish every future metadata field. |
| Native-range source and terrestrial/projected masks | Geographic projection domain | The exact source polygon and processing rules, including any Map of Life input used, require confirmation. |
| 32 controlled and 32 raw GF batch objects, cumulative curves and RDAforest fit/spatial objects | Full-panel nonlinear models and projections | The objects are large and remain external. A public manifest must identify their checksums, model settings and access route. |
| `18_LFMM2_group_locality_scores_v2.csv` and `32_v4_all_group_locality_ALT_frequency.csv` | Exact Figure 2 v5 regeneration at 185 sampled localities | These contain precise locality information. Deposition requires a sensitivity and permissions review or a controlled access route. The figure and aggregated decile summary alone do not reproduce every panel. |
| *Ameiva* assembly GCA_046270845.1 and *Aspidoscelis tigris* GCA_023333525.2 genome/GFF | Liftoff transfer and SNP annotation | The assembly accessions are recorded; Downloaded file versions, checksums and redistribution terms require documentation. |
| Manuscript text and supplementary drafts | Editorial work | The unpublished manuscript is not distributed here. |

The present repository includes selected small group, effect and diagnostic tables, plus a display-only country-border file. The tracked tables do not include an individual genotype matrix or a precise sampled-locality table. A RADseq data accession has not been verified in the available records.
