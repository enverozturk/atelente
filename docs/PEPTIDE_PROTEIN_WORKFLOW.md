# Peptide/protein workflows

These commands analyse peptide quantities relative to protein or protein-group quantities. The neo/spanning workflow is described in [Usage](USAGE.md).

## Commands

Run from the source-checkout root with a configuration file:

| Command | Purpose | Configuration example |
| --- | --- | --- |
| `Rscript scripts/normalize.R --config normalize.json` | Sample peptide/protein ratios | [MaxQuant](../inst/extdata/config/maxquant.example.json), or the generic specification below |
| `Rscript scripts/quantify.R --config quantify.json` | Descriptive protein-adjusted effects from a mapping run | [Quantification](../inst/extdata/config/quantification.example.json) |
| `Rscript scripts/statistics.R --config statistics.json` | Tests of exported sample ratios | [Statistics](../inst/extdata/config/statistics.example.json) |
| `Rscript scripts/refresh_report.R --config report.json` | Attach quantities/tests to a mapping report | [Report](../inst/extdata/config/report.example.json) |

Paths resolve relative to the configuration. Use a new output directory. Tables can be CSV, TSV/TXT or XLSX. CSV accepts comma or semicolon delimiters. Numeric text uses decimal points without thousands separators. Excel support requires readxl/openxlsx; statistical tests require limma.

## Sample ratios from linked tables

Supply peptide, protein and sample tables, with exact source column names:

```json
{
  "peptides": "peptides.csv",
  "proteins": "proteins.csv",
  "samples": "samples.csv",
  "output_dir": "ratios_results",
  "excel": true,
  "spec": {
    "observation_level": "peptide_observation",
    "link_separator": ";",
    "peptide": {
      "id": "Observation",
      "sequence": "Sequence",
      "protein_links": "Protein groups",
      "modifications": "Modification",
      "scale": "linear",
      "quantity_type": "Reported peptide intensity"
    },
    "protein": {
      "id": "Group",
      "accessions": "Members",
      "entity_type": "protein_group",
      "scale": "linear",
      "quantity_type": "Reported protein-group abundance"
    }
  }
}
```

Sequence/links, protein ID, quantity scales and quantity-type descriptions are required. Peptide ID, modification and contaminant mappings are optional, as are protein accession and contaminant mappings. Missing peptide IDs are generated for the import; original rows and modification variants remain distinct. Protein IDs must be nonempty and unique.

`entity_type` is `protein` or `protein_group`. `link_separator` separates candidate entity IDs; use `null` for one literal ID. Link whitespace is preserved unless `trim_link_whitespace` is true. Repeated links are deduplicated, and unresolved links remain diagnostic.

The Samples table uses these exact columns:

| Sample | Condition | Use | Peptide column | Protein column |
| --- | --- | --- | --- | --- |
| Control_1 | Control | Yes | Intensity 1 | Abundance 1 |
| Treatment_1 | Treatment | Yes | Intensity 2 | Abundance 2 |

Sample names and included quantity columns must be unique. `Use` is `Yes` or `No`; conditions can be blank or identical for descriptive ratios. Select worksheets with `peptide_sheet`, `protein_sheet` and `sample_sheet`.

Each layer declares `scale` as `linear` or `log2`. The sample ratio is `peptide_log2 - protein_log2`, using positive finite linear values or finite log2 values. Measurements are used without imputation or pseudocounts. Linear ratio overflow/underflow is flagged while a finite log2 ratio remains available.

Contaminant flags exclude affected observations and retain their reasons. A contaminated linked entity excludes its observation. Multiple candidate entities provide alternative denominators for the same observation.

Outputs include `sample_ratios.csv`, source/canonical tables, `normalization_manifest.json` and optional `sample_ratios.xlsx`. The R interface is `pm_import_quant_tables(peptides, proteins, spec, samples)` followed by `pm_sample_ratios(data)`.

## MaxQuant input

Use `adapter: "maxquant_peptides"` for sequence-summary `peptides.txt` and `proteinGroups.txt` from the same analysis. The [example configuration](../inst/extdata/config/maxquant.example.json) specifies the preset. Select matching experiment suffixes in the Samples table.

Peptide `Protein group IDs` link to protein-group `id`. Peptide quantities use experiment-level Intensity. Protein quantities default to Intensity, with explicit LFQ intensity or iBAQ options; record the chosen abundance measure. Group member accessions are metadata.

Contaminant, reverse and site-only flags are recorded separately. A protein group identified only by site is excluded as a denominator by default. The input observation level is sequence-summary, with unavailable modification fields retained as unavailable.

## Descriptive effects from a mapping run

`scripts/quantify.R` uses `results_dir`, `samples` and `output_dir`. Mapping artifacts may be in a flat result folder or its `audit/` subfolder. Record the protein abundance source in the configuration.

The sample workbook contains Samples and Comparisons sheets. Samples uses the five columns shown above; Comparisons uses `Name`, `Treatment` and `Control`.

For each mapped master accession, an exact protein-table accession supplies a candidate baseline. Duplicate rows, missing baselines and contaminants retain unavailable states. Multiple baselines remain alternatives.

Within each condition, both peptide and protein must have positive finite quantities. The adjusted log2 effect is the peptide condition change minus the protein condition change on those common samples. This descriptive command exports `protein_adjusted_effects.xlsx` with quantities, effects and settings.

## Statistical tests

`scripts/statistics.R` requires `ratios_dir`, `samples`, `comparisons`, `output_dir` and `design: "independent_units"`. It checks the recorded normalization and ratio hashes before use.

Sample metadata uses `Sample`, `Condition`, `Use` and `Experimental unit`. Sample IDs refer to the ratio export and selected experimental units must be unique. Existing nonblank conditions must agree with the new metadata. Comparisons uses `Name`, `Treatment` and `Control`. XLSX defaults to Samples and Comparisons sheets; `sample_sheet` and `comparison_sheet` override those names.

The response is the log2 peptide/protein ratio. An independent-sample OLS model with limma moderation estimates the treatment-minus-control difference. Tests use observations with one declared protein/group candidate and require at least two usable units per condition and two eligible observations for moderation. Insufficient coverage, exact zero residual variance, ambiguous denominators and numerical failures retain diagnostic statuses.

BH adjustment uses finite tested observation p-values within each contrast. The fit uses observed values without imputation. Output includes `results.csv`, `measurements.csv`, sample/comparison tables, `statistics_manifest.json` and optional `peptide_statistics.xlsx`. Infinite prior degrees of freedom are recorded by `prior_df_infinite`. The R interface is `pm_test_ratios(ratios, samples, comparisons, design = "independent_units")`.

A synthetic example is available:

```sh
Rscript examples/statistics_demo/generate.R my_example
Rscript scripts/normalize.R --config my_example/normalize.json
Rscript scripts/statistics.R --config my_example/statistics.json
```

## Attach quantities and export a report

`scripts/refresh_report.R` takes the original mapping report and a normalization export from the same complete, ordered source peptide table. Supply `mapping_dir`, `ratios_dir`, optional `statistics_dir` and `output_dir`. The mapping report must have an empty quantity layer; optional tests must come from that normalization export.

The attachment checks source-table order, row IDs, sequences, modifications and recorded hashes. It accepts the mapping export's `peptide_id`/`source_row` prefix after identity checks and retains the existing mappings, ratios and tests. Keep source HTML and audit artifacts together.

To regenerate a compact result bundle:

```sh
Rscript tools/refresh_output_bundle.R INPUT_REPORT NEW_DIRECTORY
```

The bundle includes the existing measurements and fitted results. Output is staged in a temporary sibling directory and moved to the destination on success. Large HTML conversions may require more memory than a fresh run; rerunning the original inputs is another way to produce a current-format bundle.
