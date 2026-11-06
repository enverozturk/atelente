# Using Atelente

Start with the [README](../README.md) to install and open the local Shiny interface. This guide covers input preparation, downloads and scripted use. [Analysis rules](ANALYSIS_RULES.md) define the reported quantities and tests.

## Prepare the input

Use a peptide table and the FASTA associated with its search. Atelente maps reported sequences to the supplied FASTA rather than searching raw spectra. Different database releases can change mappings and coordinates.

The Shiny workflow accepts CSV, TSV/TXT and XLSX. Select a worksheet when needed. Headers must be unique and on the first row. Select the sequence, master/representative accession set and sample quantity columns. Add modification and contaminant columns when available. A missing modification column is not treated as evidence of an unmodified form.

Select linear or log2 quantities explicitly. Use decimal points in numeric text, with no thousands separators. Missing measurements can be left empty. The program preserves zeros separately and uses only positive finite linear quantities for ratios. Malformed selected quantities stop the run with row/column information.

Each selected column should represent the sample whose name and condition you enter. For statistical tests it must be an independent experimental unit. Technical injections or condition averages are not additional independent samples. Atelente does not infer paired designs or merge replicates. A comparison is optional. Without one, quantities remain available and comparison effects/tests are unavailable.

The standard Shiny profile is trypsin. Select lysine blocking only when it applies to the preparation. If a pre-digestion N-terminal label applies, enter its modification name as represented in your input. These settings affect candidate classification. They are not inferred from the file name or quantity headings.

Contaminant information is optional. When supplied, flags, listed accessions and contaminant-sequence matches are used for exclusions and remain traceable in the audit. Absence of such information does not trigger a requirement to supply a list.

## Read and move a run

The default location is `atelente_runs/<run>/results/` under the R working directory. The setup allows a different run location. Keep the generated folder together so report downloads can reach adjacent files.

- **Results:** representative neo/spanning comparisons and available tests.
- **Peptides:** candidate and supporting forms, modifications and sample quantities.
- **Settings:** input mappings, samples, comparison and run settings.

In the HTML, **Download current view** exports the filtered rows across all pages. **Choose tables** provides optional All pairs, Other candidates and Annotations worksheets. Filtered exports record the current display settings and preserve the original calculated effects and tests. **Download full audit** retrieves the separate audit ZIP.

The sample counts mean measured/total within each condition. They do not pair a control replicate with a treatment replicate. Selecting a result opens supporting forms and the protein-position map. Alternative pairs are available without changing the fixed representative used in the overview.

## Use the installed package from R

The five example configurations are installed with the package:

```r
example <- system.file("extdata", "config", "config.example.json",
                       package = "atelente", mustWork = TRUE)
file.copy(example, "my_config.json", overwrite = FALSE)
```

Edit the copied file before running it. Its paths are placeholders. The basic example describes mapping and does not itself configure a full neo/spanning comparison. Remove optional protein/worksheet entries if you do not supply them. Choose an empty output directory.

After editing:

```r
config_path <- normalizePath("my_config.json", mustWork = TRUE)
config <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
atelente::pm_run(config, config_base = dirname(config_path))
```

Paths resolve relative to the configuration file. Use forward slashes in Windows JSON paths. The Shiny run's saved configuration records its explicit column/sample setup and can be used to reproduce that setup with the associated input files.

## Use a source checkout

From the repository root, an edited configuration can be run with:

```sh
Rscript run.R --config my_config.json
```

The [mapping configuration example](../inst/extdata/config/config.example.json) is also in the source tree. Scripts in the repository root are source-checkout tools, not files installed on the executable search path by the R package.

A small existing example illustrates linked mapping and the separate peptide/protein quantity workflow:

```sh
Rscript examples/linked_report/generate.R atelente_demo
```

It requires limma and a new output directory, uses synthetic data and makes no reference downloads. Open `atelente_demo/report/report.html`. This example produces a linked report rather than the full current Shiny output bundle.

## MaxQuant peptide/protein quantities

The separate `maxquant_peptides` adapter reads sequence-summary `peptides.txt` and `proteinGroups.txt` from the same analysis. Use the [MaxQuant configuration example](../inst/extdata/config/maxquant.example.json) and a Samples CSV/XLSX with `Sample`, `Condition`, `Use`, `Peptide column` and `Protein column`. Select matching experiment suffixes explicitly, for example `Intensity Sample1` in both tables.

```sh
Rscript normalize.R --config my_maxquant.json
```

This command calculates peptide/protein-group ratios. It does not run the neo/spanning workflow. The preset links literal `Protein group IDs` to the protein table's `id`. Peptide quantities are experiment-level Intensity. The protein denominator defaults to Intensity, with explicit LFQ intensity or iBAQ options. These quantity types are not interchangeable estimators.

Sequence-summary rows are not modification-specific measurements. Missing modification information stays unavailable. Optional contaminant, reverse and site-only flags have separate meanings. A group with `Only identified by site` is excluded as a denominator by default. This preset is not a claim of universal MaxQuant-version or DIA-format support.
