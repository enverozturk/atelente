# Mapping-dependent protein adjustment

`quantify.R` calculates descriptive peptide/protein-adjusted condition effects from an existing mapping run. The main neo/spanning analysis does not require this separate command.

Supply `results_dir`, `samples` and a new `output_dir` using the installed [quantification configuration](../inst/extdata/config/quantification.example.json). Record the source of the protein abundances. Required mapping artifacts are resolved from either the historical flat output or the current `audit/` folder.

```sh
Rscript quantify.R --config my_quantification.json
```

The sample workbook needs Samples and Comparisons sheets. Samples uses `Sample`, `Condition`, `Use`, `Peptide column` and `Protein column`. Comparisons uses `Name`, `Treatment` and `Control`. Select existing columns explicitly. Raw, normalized and grouped representations are distinct and are not substituted silently.

For each mapped master accession, an exact protein-table accession supplies a possible baseline. Duplicate protein rows, missing baselines and contaminants retain unavailable states. Multiple baselines are alternatives for one observation.

Within each condition, calculations use samples with both peptide and protein quantities positive and finite. The adjusted log2 effect is the peptide condition change minus the protein condition change on those same samples. No imputation, replicate pairing or p-values are introduced by this command.

The output includes `protein_adjusted_effects.xlsx` with traceable quantities, effects and settings. Input paths and original source identities remain available in the run artifacts. The generic alternative that starts from explicitly linked tables is documented in [Generic input](GENERIC_INPUT.md).
