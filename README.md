# Atelente

**Analyse neo-terminal peptides relative to peptides spanning the same cleavage bond.**

Atelente is an R package with a local Shiny interface for shotgun proteomics. It maps observed peptides to protein sequences, links neo-N, neo-C and spanning forms, and compares within-sample neo/spanning ratios between conditions. Its primary application is data collected without terminal-peptide enrichment.

## Installation

Use R 4.3 or later:

```r
install.packages(c("remotes", "shiny", "readxl", "openxlsx"))
remotes::install_github("enverozturk/atelente", dependencies = NA)
```

For statistical testing, install limma:

```r
install.packages("BiocManager")
BiocManager::install("limma", ask = FALSE, update = FALSE)
```

The repository currently requires an account with access for installation.

## Quick start

```r
shiny::runApp(atelente::pm_wizard_app())
```

1. Choose a peptide table and its matching FASTA.
2. Select the sequence, protein accession and sample quantity columns.
3. Enter sample names and conditions; select the comparison and digestion settings.
4. Review the settings and run the analysis.

Use decimal points in numeric text and leave missing measurements empty. For statistical testing, each sample must be an independent experimental unit. Prepare technical-replicate summaries before import.

MEROPS and UniProt annotations are optional. Their first use downloads reference records; later runs reuse the local cache. Input tables are processed on your computer.

## Results

| File | Contents |
| --- | --- |
| `report.html` | Interactive results, supporting peptides, protein maps and reference annotations |
| `atelente_results.xlsx` | Results, Peptides and Settings worksheets |
| `atelente_audit.zip` | Measurements, alternative pairs, mapping and provenance |
| `run_manifest.json` | Run settings, reference details and file hashes |

Results are saved under `atelente_runs/<run>/results/` by default. Keep the result files together and open the HTML report in Chrome or Edge.

The reported effect is the treatment-minus-control difference in mean log2(neo/spanning). Positive values indicate a higher ratio in treatment. Optional limma tests use independent samples, with Benjamini–Hochberg adjustment within each comparison and test family.

## Documentation

- [Usage](docs/USAGE.md): input preparation, Shiny and scripted runs.
- [Analysis](docs/ANALYSIS_RULES.md): candidate bonds, pair selection and statistics.
- [Reference annotations](docs/REFERENCE_ANNOTATIONS.md): MEROPS, UniProt and offline use.
- [Peptide/protein workflows](docs/PEPTIDE_PROTEIN_WORKFLOW.md): additional normalization, statistics and report commands.

Source-checkout commands are in `scripts/`, runnable examples in `examples/`, and bundled example tables/configurations in `inst/extdata/`. Package functions, help pages and tests are in `R/`, `man/` and `tests/`.

## Study data and citation

[Study inputs and analysis records](https://drive.google.com/file/d/1a0HeQokmYq-FII6ryM4Ctrw24Z6FDpXq/view) are archived separately; the working archive currently requires access. The manuscript analyses used [commit b20990f](https://github.com/enverozturk/atelente/tree/b20990f7243c87fa468220c7ad27a5af2ce0559a). Software citation details are in [CITATION.cff](CITATION.cff).

GPL-3. Author and maintainer: Enver Ozturk ([ORCID](https://orcid.org/0009-0004-8359-8032)); [enver10166@gmail.com](mailto:enver10166@gmail.com).
