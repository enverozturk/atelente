# Usage

Install Atelente using the [README](../README.md). The Shiny app provides four steps: Files, Columns, Samples, and Review & run.

## Prepare the input

Supply a peptide table and the FASTA used for its database search. Coordinates refer to that FASTA. Accepted table formats are CSV, TSV/TXT and XLSX; select a worksheet for Excel input. Put unique column headers on the first row.

Select the peptide sequence, master/representative accession set and sample quantity columns. Modification and contaminant columns are optional. An unavailable modification state remains recorded as unavailable.

Declare whether quantities are linear or log2. Use decimal points without thousands separators. Leave missing measurements empty. Ratios use positive finite linear quantities; zeros and missing values are recorded separately. An invalid value in a selected quantity column produces a message identifying its row and column.

Give each included quantity column a sample name and condition. Statistical comparisons use independent experimental units. Summarize technical replicates before import. A run can also retain measurements without a treatment/control comparison.

## Set up and run

The standard Shiny digestion profile is trypsin. Select lysine blocking when applicable and enter any pre-digestion N-terminal label using its name in the input. These settings determine candidate classification.

Contaminant annotations can come from input flags, accession lists or sequence matches. Excluded observations and their reasons remain in the audit.

Reference annotations are enabled by default. The first run downloads the required records unless a cache or supplied snapshots are available. See [Reference annotations](REFERENCE_ANNOTATIONS.md) for offline use.

Choose a destination and run the analysis. The default is `atelente_runs/<run>/results/` under the R working directory.

## Read the results

The Excel workbook has three default sheets:

- **Results:** representative neo/spanning comparisons and tests. `Spanning ends` and `Spanning terminal form` show whether the spanning form is itself a terminal form (see [Analysis](ANALYSIS_RULES.md#representative-pair-selection)).
- **Peptides:** observed forms, modifications and sample quantities.
- **Settings:** selected columns, samples, comparison and run settings.

In the HTML report, select a result to inspect its peptides and protein map. Sample counts show measured/total observations within each condition. Alternative pairs are available alongside the fixed representative.

**Download current view** exports the filtered rows across all pages. **Choose tables** adds optional All pairs, Other candidates and Annotations worksheets. **Download full audit** retrieves the audit ZIP. Display filters preserve the fitted effects, tests and representative selection.

Keep the result folder together. Extract downloaded ZIPs before opening the report so its adjacent downloads remain accessible.

## Use the installed package from R

Copy an example configuration:

```r
example <- system.file("extdata", "config", "config.example.json",
                       package = "atelente", mustWork = TRUE)
file.copy(example, "my_config.json", overwrite = FALSE)
```

Edit its input, output and column settings. The basic configuration runs peptide mapping. A saved Shiny configuration also contains the selected sample, digestion and neo/spanning settings. Use a new output directory.

```r
config_path <- normalizePath("my_config.json", mustWork = TRUE)
config <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
atelente::pm_run(config, config_base = dirname(config_path))
```

Paths resolve relative to the configuration file. Use forward slashes in Windows JSON paths. Optional protein-table and worksheet fields can be removed when unused.

## Use a source checkout

From the repository root:

```sh
Rscript scripts/run.R --config my_config.json
```

The [configuration examples](../inst/extdata/config) are included in the checkout. The root `app.R` also supports `shiny::runApp(".")` from the checkout.

For a runnable example of mapping linked to peptide/protein quantities:

```sh
Rscript examples/linked_report/generate.R atelente_demo
```

This example uses synthetic data, requires limma and writes `atelente_demo/report/report.html`. It illustrates the [peptide/protein workflow](PEPTIDE_PROTEIN_WORKFLOW.md).
