# RNA-seq pipeline — GMCF project 3818

MC3T3-E1 pre-osteoblasts, vehicle vs 300 nM SB273005.
Quantification with kallisto, differential expression with DESeq2.

## Design

| Group | n | Libraries |
|---|---|---|
| Vehicle | 5 | Samples 1, 3, 4, 5, 8 |
| SB273005 300 nM | 6 | Samples 9, 12, 13 + three replacements |
| No-template control | 1 | excluded from statistics |

Eighteen RNA samples came from a single lysis and a single Qiagen RNeasy
extraction. Twelve were submitted; ten were selected for sequencing. Samples 10
and 11 failed library preparation at the core across five attempts, producing
libraries the core described as resembling negative controls. Their RIN (10) and
purity ratios were comparable to samples that prepped successfully, so no cause
was identifiable from the QC data. Three replacements from the same RNA
extraction batch were submitted and prepped later. This is the source of the
library-prep batch variable in the model.

## Pipeline

| Script | Step |
|---|---|
| `01_strandtest.sh` | Determine strandedness empirically |
| `02_runall.sh` | Quantify all libraries with kallisto |
| `03_makecounts.py` | Aggregate transcripts to genes |
| `04_makefigs.py` | Normalisation, PCA, gene-level figures |

Requires kallisto >= 0.50 and the prebuilt mouse transcriptome index from
pachterlab/kallisto-transcriptome-indices.

## Strandedness

Determined from the data rather than assumed, by quantifying one library three
ways and comparing pseudoalignment rates:

| Setting | Pseudoaligned |
|---|---|
| unstranded | 92.3% |
| --fr-stranded | 1.0% |
| --rf-stranded | 88.8% |

Unstranded scores highest because it accepts both orientations; the informative
comparison is forward vs reverse. The libraries are reverse-stranded (dUTP
protocol), so `--rf-stranded` was used throughout.

## Normalisation

Two tables are produced and they are not interchangeable.

- **Raw counts** are the input to DESeq2. DESeq2 models counts with a negative
  binomial distribution and computes its own size factors (median-of-ratios)
  internally, so pre-normalised input would invalidate the model.
- **TPM** is corrected for library size and transcript length, and is used for
  reading off expression levels.

Library sizes before normalisation ranged 28.0-40.2 M (spread 8.9%); after
size-factor correction, 31.9-34.1 M (spread 1.9%).

## Model

    ~ library_batch + treatment

Batch is included because principal component analysis showed the two library
preparations separating completely on PC1 (49% of total variance), with
treatment groups intermixed. A direct batch contrast found 7,147 of 13,138 genes
(54%) differing at FDR 5%, including Actb at 1.6-fold — technical, not
biological. Treatment and batch are not confounded, because the original
preparation contains both treatment groups; the comparison was also run on that
batch alone as a cross-check.

Genes were retained with >= 10 counts in >= 3 samples, leaving 13,138 tested.

## Result

No gene reached significance at FDR 5%. Smallest adjusted p-value 0.14.

The null is well supported rather than underpowered:

- median standard error on log2 fold change 0.071 across 8,026 well-expressed
  genes, giving ~80% power to detect a 1.15-fold change
- median observed fold change 1.03x; largest 1.80x; none above 2x
- 78 genes at raw p < 0.01 against 131 expected by chance

Interpretation: SB273005 alters migration without altering transcription, so the
effect is post-translational on integrin signalling.

Incidental finding: Itgb5 is expressed 33-fold higher than Itgb3 in these cells
(573 vs 17 TPM), so the drug acts predominantly through alpha-v beta-5.

## Library QC

2 x 150 bp paired-end, single lane, unique dual indexes. 33.8-43.3 M read pairs
per library, 93-96% of bases at Q30, pseudoalignment 82.7-92.9%, no adapter
contamination, no flowcell tile below 90% Q30, RIN 9.8-10. The no-template
control reached 1.9% pseudoalignment, confirming negligible cross-contamination.

## Limitations

- Single timepoint; a transient response would be invisible
- Single dose sequenced (300 nM); the 1 uM arm was not sequenced
- With batch in the model the treatment estimate leans on the original batch
- No protein-level data; phospho-westerns are the appropriate follow-up
