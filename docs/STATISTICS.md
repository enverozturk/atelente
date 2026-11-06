# Separate peptide/protein statistics

This command tests peptide/protein log2 ratios exported by [normalize.R](GENERIC_INPUT.md). It is separate from the main [neo/spanning pair model](PAIR_STATISTICS.md).

## Input contract

Provide `ratios_dir`, `samples`, `comparisons`, `output_dir` and `design: "independent_units"`, following the installed [statistics configuration](../inst/extdata/config/statistics.example.json). Paths resolve relative to that configuration. The recorded normalization and ratio hashes are checked before use.

The sample metadata has exactly `Sample`, `Condition`, `Use` and `Experimental unit`. Sample IDs refer to the ratio export. Selected experimental units must be unique. An existing nonblank condition must agree with the new metadata. The comparison table has exactly `Name`, `Treatment` and `Control`.

```sh
Rscript statistics.R --config my_statistics.json
```

Use a new output directory. XLSX inputs default to Samples and Comparisons worksheets. `sample_sheet` and `comparison_sheet` can override them. limma is required. Excel input/output also needs readxl/openxlsx.

## Model and interpretation

The response is a sample log2 peptide/protein ratio. The effect is the treatment mean minus the control mean. An independent-sample OLS model with limma moderation is fitted to observations with one declared protein/group candidate. Alternative denominators do not create independent tests or establish unique protein attribution.

This separate model requires at least two usable units per condition and two eligible observations for moderation. Exact zero residual variance, insufficient values, ambiguity and numerical failures remain diagnostic. Missing values are not imputed. Paired, repeated, batch and factorial models are not fitted. Technical repeats cannot create independent units by relabelling them.

BH adjustment is within each contrast's finite tested observation p-values. Filters do not redefine this family. Minimum counts do not guarantee power or stable moderation. Condition-dependent missingness and shared protein denominators can affect inference.

## Outputs

`results.csv` retains descriptive effects, test results and diagnostics. `measurements.csv` records sample values and fit use. Samples/comparisons, `statistics_manifest.json` and optionally `peptide_statistics.xlsx` retain design and provenance. Infinite prior degrees of freedom use an explicit `prior_df_infinite` flag.

The R API is `pm_test_ratios(ratios, samples, comparisons, design = "independent_units")`. Existing results can be attached to a matching report using [refresh_report.R](REPORT_QUANTITIES.md) without refitting.

A synthetic developer example is available from a source checkout:

```sh
Rscript examples/statistics_demo/generate.R my_example
Rscript normalize.R --config my_example/normalize.json
Rscript statistics.R --config my_example/statistics.json
```

This is a software example. Its results are not an experimental benchmark or a general calibration claim.
