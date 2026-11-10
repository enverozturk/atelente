# Mouse kidney — MSV000094661

[TermineR study](https://doi.org/10.1002/pmic.202300491)

Run from the repository root:

```sh
Rscript manuscript/run.R 02_Kidney_MSV000094661
```

The input contains 28,710 prepared peptidoforms, with six WT and five KO quantity columns. The saved settings use trypsin specificity with lysine blocking, the TMTpro pre-digestion N-terminal label, representative/master accession assignment and selection k = 3. The comparison is KO minus WT, with limma tests and Benjamini–Hochberg correction. The sample design is included directly in `analysis.json`; CSV copies and the deposited channel annotation are also supplied.

The native Results table used in the manuscript contains 6,056 rows. In the manuscript comparison, 2,704 of the 2,740 unique published terminal features were matched, 1,181 had an observed spanning peptide, and 1,125 had finite selected comparisons. Fold-change signs agreed for 847 of those comparisons (75.3%). The supplied `published_terminal_features.csv.gz` retains the source table used for that comparison; these matched comparisons are a derived subset of the native results.

The reference FASTA is byte-identical to `protein.fas` distributed with the TermineR repository (`data-raw/pkd_mouse_model_search`). Prepared quantities preserve the selected source PSM quantity vectors, without additional normalization or imputation. `metadata.json` records input hashes and source details. The manuscript analyses used limma 3.68.5.
