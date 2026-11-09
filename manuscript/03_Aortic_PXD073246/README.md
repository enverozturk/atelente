# Aortic valve — PXD073246

[Source study](https://doi.org/10.1186/s10020-026-01499-0)

Run from the repository root:

```sh
Rscript manuscript/run.R 03_Aortic_PXD073246
```

The input retains all 14,597 deposited peptide rows. Consecutive technical triplets were averaged over observed normalized linear intensities before import, giving nine Control and ten Severe samples. The study analyst confirmed this grouping. The saved settings use trypsin specificity without lysine blocking, representative/master accession assignment and selection k = 3. The comparison is Severe minus Control, with limma tests and Benjamini–Hochberg correction. No protein fold-change prefilter or imputation is applied.

The validated native Results table contains 2,165 rows, including 1,857 tested pairs and 454 with q < 0.05. Among the significant pairs, 370 have a positive ratio log2 fold change.

`sample_design.xlsx` is the design used by the configuration. CSV copies, `triplet_mapping.csv`, the deposited sample list and the community-annotated SDRF document the sample assignments. The 1,872-record reference FASTA was recovered by exact UniProt accession and is the same FASTA used in the validated Atelente analysis; it is not the original MatrixDB search FASTA. `metadata.json` records input hashes and source details. The replicated analyses were validated with limma 3.68.5.
