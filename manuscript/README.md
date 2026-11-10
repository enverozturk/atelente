# Manuscript analyses

These folders contain the prepared peptide quantities, FASTA sequences, sample information and configurations used for the three Atelente manuscript analyses. They also serve as worked examples. The original studies are linked below; raw mass-spectrometry files remain in their source repositories.

## Run an analysis

Install Atelente and its optional interface/Excel dependencies using the [main README](../README.md). Install limma for the kidney and aortic analyses. Download and extract this repository, then run one of these commands from its top-level directory:

```sh
Rscript manuscript/run.R 01_GluC_PXD000503
Rscript manuscript/run.R 02_Kidney_MSV000094661
Rscript manuscript/run.R 03_Aortic_PXD073246
```

Each command uses that dataset's `analysis.json` and writes a `results/` folder beside its inputs. The runner expands `peptides.csv.gz` to `peptides.csv`, verifies the uncompressed peptide and FASTA hashes, and calls `atelente::pm_run()`. A run stops if `results/` already contains files; move or rename it before running again.

Annotations are enabled. The first run requires internet access to retrieve MEROPS and UniProt records; subsequent runs reuse the local cache. The manuscript analyses used UniProt release 2026_04. Reference identifiers and hashes are recorded in [reference_versions.json](reference_versions.json). New retrievals may use updated reference records. [Reference annotations](../docs/REFERENCE_ANNOTATIONS.md) explains how to supply frozen snapshots for offline use.

The replicated analyses used limma 3.68.5. Use the same version to reproduce the reported statistics; the runner warns if another limma version is used. Running these inputs with Atelente v0.1.0 reproduces the Results row counts below. Each new run records its own software and reference provenance.

## Datasets

| Dataset | Prepared rows | Conditions | Results rows |
| --- | ---: | --- | ---: |
| [GluC, PXD000503](01_GluC_PXD000503) | 3,017 | 16 h versus 1 h; eight channels reported | 133 |
| [Kidney, MSV000094661](02_Kidney_MSV000094661) | 28,710 | Six WT and five KO | 6,056 |
| [Aortic valve, PXD073246](03_Aortic_PXD073246) | 14,597 | Nine Control and ten Severe | 2,165 |

Each dataset has a short README, a runnable configuration, source/preparation metadata and its sample/comparison information. [SHA256SUMS.csv](SHA256SUMS.csv) identifies the distributed files.

Results rows are the rows of the Results sheet in `atelente_results.xlsx`. The manuscript supplementary tables are derived from these outputs; for example, the 1,125 kidney comparisons are a subset of the 6,056 rows.

## Source attribution

- **GluC:** Schlage and colleagues, [doi:10.1074/mcp.M113.035139](https://doi.org/10.1074/mcp.M113.035139), dataset PXD000503. Prepared from the published supplementary quantities.
- **Kidney:** the [TermineR study](https://doi.org/10.1002/pmic.202300491), dataset MSV000094661. Prepared peptide quantities and the reference FASTA are those used in the manuscript analysis. The published terminal-feature table is supplied separately for the manuscript comparison.
- **Aortic valve:** [doi:10.1186/s10020-026-01499-0](https://doi.org/10.1186/s10020-026-01499-0), dataset PXD073246. Prepared quantities retain the source normalization; technical triplets are averaged before import.

These are attributed study inputs, separate from the Atelente software license. Cite the original studies when reusing their data. Reference FASTA files used in the analyses are included; their provenance is specified in each dataset's metadata.
