# GluC time course — PXD000503

[Source study](https://doi.org/10.1074/mcp.M113.035139)

Run from the repository root:

```sh
Rscript manuscript/run.R 01_GluC_PXD000503
```

The input retains all 3,017 observations and eight iTRAQ quantity channels from published supplementary table 1. Each time/control channel is a separate condition. The saved settings use trypsin specificity with lysine blocking, the iTRAQ pre-digestion N-terminal label, unique assignment and selection k = 1. Statistical tests are disabled.

The validated native Results table contains 133 rows. The analysis links the complementary COL5A2 products at E1185|G1186 to their spanning peptide. The 842-record reference FASTA is the one used in the validated Atelente analysis; it is not identified as the original study search database.

`sample_design.xlsx` is the sample design used by the configuration. `samples.csv` and `comparisons.csv` provide readable copies. Quantities retain zeros and missing measurements, with no additional scaling or imputation. `metadata.json` records input hashes and source details.
