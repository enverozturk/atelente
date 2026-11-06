# Atelente

**Link neo-terminal and spanning peptides at the same candidate cleavage bond.**

Atelente is an R package with a local Shiny interface for shotgun proteomics. It places observed peptides on supplied protein sequences, links neo-N, neo-C and spanning forms, and compares neo/spanning ratios between conditions when both forms are measured in the same sample. Its primary use is data without terminal-peptide enrichment.

The spanning peptide provides a reference across the same bond, complementary to protein-level abundance. Results combine searchable tables, protein-position maps and MEROPS/UniProt annotations.

> Development version. While this repository is private, GitHub installation requires authentication with an account that has repository access.

## Install and open

Use R 4.3 or later. For the Shiny interface and Excel input/output:

```r
install.packages(c("remotes", "shiny", "readxl", "openxlsx"))
remotes::install_github("enverozturk/atelente", dependencies = NA)
shiny::runApp(atelente::pm_wizard_app())
```

For optional statistical tests, also install limma:

```r
install.packages("BiocManager")
BiocManager::install("limma", ask = FALSE, update = FALSE)
```

## Analyse your data

1. **Files:** choose a peptide table (CSV, TSV/TXT or XLSX) and its matching FASTA, optionally gzip-compressed.
2. **Columns:** select the peptide sequence, master/representative protein and quantity columns. Add modifications when available and select linear or log2 quantities.
3. **Samples:** name the samples and conditions. Optionally select a treatment/control comparison and statistical tests.
4. **Review & run:** check the digestion settings and reference options, then run the analysis.

Before you start:

- Use decimal points in numbers (`100.5`).
- Leave missing measurements empty. Do not replace them with zeros.
- Tests assume independent samples. Do not enter technical replicate measurements as separate samples.

**Reference annotations.** When enabled, Atelente downloads the required MEROPS and UniProt references and saves them on your computer. Later runs work offline when those references are already available. UniProt lookups send protein accession identifiers. Your input tables are processed locally. You can turn annotations off in the setup.

## Read and download results

| Output | Contents |
| --- | --- |
| `report.html` | Searchable results, supporting peptides, alternative pairs, protein maps and reference matches |
| `atelente_results.xlsx` | Results, Peptides and Settings worksheets |
| `atelente_audit.zip` | Detailed measurements, links, exclusions and provenance |
| `run_manifest.json` | Run metadata, reference provenance and file hashes |

By default, the app saves each run in `atelente_runs/<run>/results/` inside your R working directory. The report can export the current filtered view or selected worksheets, including All pairs and Annotations.

Keep the results folder together. If you receive it as a ZIP file, extract it before opening the report in Chrome or Edge. A finished report works offline.

## What the main result means

Ratio log2FC compares the mean log2(neo/spanning) ratio between treatment and control, using samples where both forms were quantified. Positive values indicate a higher ratio in treatment. One-condition observations are reported without a numerical fold change.

A fixed rule chooses a representative pair for each bond and terminal role, independently of effect size or p-values. Alternative pairs remain available. Optional limma tests evaluate ratio differences, with Benjamini–Hochberg adjustment separately within each comparison and test family.

## Guides

[Usage and scripted workflows](docs/USAGE.md) · [Analysis rules](docs/ANALYSIS_RULES.md) · [Reference annotations](docs/REFERENCE_ANNOTATIONS.md)

## License and author

GPL-3. Author and maintainer: Enver Ozturk ([ORCID](https://orcid.org/0009-0004-8359-8032)). Contact: [enver10166@gmail.com](mailto:enver10166@gmail.com). Software citation metadata is in [CITATION.cff](CITATION.cff).
