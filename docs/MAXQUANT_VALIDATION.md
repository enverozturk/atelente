# MaxQuant quantity input

The `maxquant_peptides` adapter supports sequence-summary `peptides.txt` and `proteinGroups.txt` from the same analysis. It belongs to the separate [peptide/protein quantity workflow](GENERIC_INPUT.md). This page documents input support, not a new biological validation result.

Use the installed [MaxQuant example configuration](../inst/extdata/config/maxquant.example.json). Supply a Samples table with `Sample`, `Condition`, `Use`, `Peptide column` and `Protein column`. Select matching experiment suffixes explicitly.

```sh
Rscript normalize.R --config my_maxquant.json
```

The preset links literal peptide `Protein group IDs` to protein-group `id`. Group member accessions are metadata. Peptide quantities use experiment-level Intensity. Protein quantities default to Intensity, with explicit LFQ intensity or iBAQ options. Record which quantity method was used. These methods are not interchangeable estimators.

Optional contaminant, reverse and site-only flags are separate. A group identified only by site is excluded as a denominator by default. Sequence-summary rows are not modification-specific measurements. An absent modification field remains unavailable.

The command calculates peptide/protein ratios without FASTA placement. It does not produce neo/spanning pairs. This adapter's contract is not a claim of support for every MaxQuant version or every DIA export.
