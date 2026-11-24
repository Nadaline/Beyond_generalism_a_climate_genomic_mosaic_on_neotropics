# Liftoff transfer and SNP annotation record

## Recorded Liftoff run

The output GFF header records **Liftoff v1.6.3** and the complete command:

```text
/tmp/liftoff_venv/bin/liftoff \
  -g GCA_023333525.2_genomic.gff3 \
  -o Aspidoscelis_on_Ameiva_GCA_046270845.1.gff3 \
  -u unmapped_features.txt \
  -dir liftoff_intermediate \
  -p 4 \
  GCA_046270845.1_genomic.fna \
  GCA_023333525.2_genomic.fna
```

The reference annotation and genome were from *Aspidoscelis tigris* assembly
GCA_023333525.2; the target genome was the *Ameiva ameiva* assembly
GCA_046270845.1 used to align the RADseq VCF.

The only computational option that differed from Liftoff's defaults was
`-p 4` (four parallel processes; default is one). The command did **not** use
`-a`, `-s`, `-copies`, `-exclude_partial`, `-polish`, or custom
`-mm2_options`.

## Transfer yield and thresholds

| Quantity | Value |
|---|---:|
| Source gene features | 19,975 |
| Gene features written to target GFF | 18,157 |
| IDs in `unmapped_features.txt` | 1,818 |
| Source genes represented in target GFF | 90.90% |
| Target-GFF genes with coverage >= 0.50 and sequence identity >= 0.50 | 9,509 |
| Proportion of source genes passing both default thresholds | 47.60% |
| Proportion of target-GFF genes passing both default thresholds | 52.37% |

Because `-exclude_partial` was not supplied, the target GFF retains mappings
below either default threshold and marks them with `partial_mapping=True` and
or `low_identity=True`. Thus, 18,157/19,975 is the **annotation-transfer
yield**, whereas 9,509/19,975 is the subset satisfying both default mapping
criteria. The target GFF contains 7,950 `partial_mapping=True` gene features
and 8,648 `low_identity=True` gene features; these categories can overlap.

Liftoff v1.6.3 defaults are a minimum child-feature alignment coverage of
`-a 0.50` and a minimum child-feature sequence identity of `-s 0.50`.
Gene-level `coverage` and `sequence_ID` values are retained in the ninth GFF
column.

## Extra gene copies

`-copies` was not supplied. Liftoff therefore did not search the target genome
for additional unannotated copies after the primary transfer. All 18,157
target-GFF gene features have `extra_copy_number=0`. Multiple mRNAs under a
gene are transferred reference isoforms, not separately inferred extra gene
copies.

## SNP-to-annotation rule

SNPs were represented as one-base genomic ranges and compared with **gene**,
exon and CDS features using `findOverlaps(..., ignore.strand = TRUE)`.

1. A SNP overlapping any CDS was labelled `CDS`; otherwise any exon was
   labelled `exonic_non_CDS`; otherwise any gene was labelled `intronic`;
   otherwise it was `intergenic`.
2. Gene identity was assigned at the **parent-gene level**, not at the
   transcript level. Consequently, a SNP overlapping multiple isoforms of the
   same parent gene received one gene-level annotation; no transcript ID was
   retained.
3. For an intergenic SNP, `distanceToNearest(..., ignore.strand = TRUE)`
   assigned the nearest gene and recorded its distance in base pairs.
4. The original script did not implement an explicit biological tie-break for
   a SNP overlapping two distinct gene features. Its vector assignment
   `gene_for_locus[queryHits(gene_hits)] <- subjectHits(gene_hits)` retains the
   last overlapping gene in imported GFF order. This is reproducible from the
   script but should be described as an ordering rule, not as evidence that
   one transcript or gene is preferentially affected.

## Sources inspected

- `work/Aspidoscelis_on_Ameiva_GCA_046270845.1.gff3` (header, gene attributes,
  transfer counts)
- `work/unmapped_features.txt`
- `work/GCA_023333525.2_genomic.gff3`
- `scripts/prepare_liftoff_ameiva.sh`
- `RDA-2027-bonferroni.R`, annotation block around lines 1622--1685

Liftoff default thresholds and option semantics are documented in the
[official Liftoff README](https://github.com/agshumate/Liftoff/blob/master/README.md).
