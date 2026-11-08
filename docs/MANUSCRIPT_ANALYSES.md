# Manuscript analyses

The manuscript analyses used Atelente commit [`b20990f7243c87fa468220c7ad27a5af2ce0559a`](https://github.com/enverozturk/atelente/tree/b20990f7243c87fa468220c7ad27a5af2ce0559a). Version 0.1.0 is a release candidate: its runtime files in `R/` and `inst/` are unchanged from that commit. The version update does not change or rerun the reported analyses. Cite the analysed commit when referring to the manuscript results.

**Archive status:** the Drive links below identify the current working archives and may require permission. They are not permanent public deposits. Public archive links, the release date and a persistent release identifier remain to be completed before publication. Archived wrappers preserve the original execution record and can require local path adaptation; this document does not claim they execute unchanged after downloading.

## Data and recorded analyses

All three native analyses used `cleavage.analysis_scope = "specificity_v3"`. Full column mappings, reference settings and sample designs are retained in the configurations below. Paths are relative to the root of the named ZIP, not this repository.

| Study and source | Design used in the manuscript | Working archive |
| --- | --- | --- |
| Schlage et al., 2014, *Time-resolved Analysis of the Matrix Metalloproteinase 10 Substrate Degradome*. [DOI: 10.1074/mcp.M113.035139](https://doi.org/10.1074/mcp.M113.035139); [PXD000503](https://www.ebi.ac.uk/pride/archive/projects/PXD000503). | The publication's **GluC** Supplementary Table S1 time course; pre-digestion iTRAQ, lysine-blocked digestion; one observation per time point, descriptive 16 h versus 1 h comparison, testing disabled. The source article title does not identify the treatment of this selected experiment. | [Schlage folder](https://drive.google.com/drive/folders/1tUZF4n0LTOOFr77VFVY_JyqOId3nLQAR): `Schlage_GluC_analyses_and_results_b20990f.zip` |
| Cosenza-Contreras et al., 2024, *TermineR: Extracting information on endogenous proteolytic processing from shotgun proteomics data*. [DOI: 10.1002/pmic.202300491](https://doi.org/10.1002/pmic.202300491); MassIVE accession **MSV000094661**. | Mouse kidney, six WT and five KO samples; pre-digestion TMTpro, lysine-blocked digestion; KO versus WT. The prepared input retains complete observed PSM quantity vectors, with no additional normalization or imputation. | [PKD folder](https://drive.google.com/drive/folders/1pPZxu7pujPi9Z0fEqMapiBNJkzho_F2g): `PKD_analyses_and_results_b20990f.zip` |
| Iacobucci et al., 2026, *Proteomic and degradomic signatures of extracellular matrix remodeling in calcific aortic valve stenosis*. [DOI: 10.1186/s10020-026-01499-0](https://doi.org/10.1186/s10020-026-01499-0); [PXD073246](https://www.ebi.ac.uk/pride/archive/projects/PXD073246). | Non-enriched semi-tryptic **DDA** peptide table; nine Control and ten SEVERE QC-supported triplet means, not named donors. Finite positive quantities were averaged before native ratio testing. Inference is conditional on the F-column grouping; columns were matched by exact header. | [Aortic folder](https://drive.google.com/drive/folders/1nX10MUxSrMZgdSkhxwhSEgOYwDtxc6CG): `Aortic_analyses_and_results_b20990f.zip` |

### Schlage GluC

- Recorded run configuration: `Schlage_GluC/config_used.json`; existing relative-path configuration: `Schlage_GluC/config_portable.json`.
- Inputs: `Schlage_GluC/inputs/Schlage_peptides_from_audit.csv` and `Schlage_GluC/inputs/Schlage_reference.fasta`.
- Original run wrapper: `Schlage_GluC/provenance/run_schlage_frozen_b20990f.R`.
- Native outputs and run manifest: `Schlage_GluC/results/`; reference snapshots and their provenance: `Schlage_GluC/references/`.

The recovered peptide table is a re-export of archived observations, not a claim of byte identity with the publication's original input file. The configuration and validation records document that distinction.

### PKD and the published-effect/reference comparisons

- Recorded native configuration: `PKD/config_used.json`; input base: `PKD/provenance/input_snapshot/`.
- Original wrapper: `PKD/provenance/run_pkd_frozen_b20990f.R`; native outputs: `PKD/results/`; execution notes: `PKD/RUN_README.txt`.
- Published S2 comparison: `PKD_comparisons/reproducibility/compare_s2_effects.R`.
- Controlled reference comparison: `PKD_comparisons/reproducibility/compare_pkd_references.R`; model and cohort record: `PKD_comparisons/controlled/analysis_manifest.json`.

The first comparison retains the reported S2 effects. The second compares an observed-peptidoform sum (proxy) with a spanning reference on a common feature/sample set and model; it is not an exact reproduction of TermineR's normalization. Original test families and q values are retained in the archived outputs.

### Aortic valve and technical-replicate preparation

- Recorded configuration and native wrapper: `Aortic/config_used.json`, `Aortic/run_analysis.R`.
- Inputs: `Aortic/inputs/subject_means.csv`, `subject_design.xlsx`, `technical_measurement_counts.csv`, `triplet_mapping.csv` and `recovered_reference_20261104.fasta` (all under `Aortic/inputs/`).
- Native outputs: `Aortic/results/`; extracted audit: `Aortic/audit/`; numerical checks: `Aortic/validation/`.
- The pre-package preparation scripts are currently in **`Manuscript_extra_scripts_REVIEW_BEFORE_UPLOAD.zip`**, in the [additional-script working folder](https://drive.google.com/drive/folders/1VuyxBtXz25KRWUJ_QEBxyicfRARWLWSd): `Scripts/aort_triplet_qc.R`, `Scripts/prepare_aort_qc_positive.R` and `Scripts/prepare_aort_inference.R`. They are not presently included in the aortic ZIP. This is an archived provenance copy, not a separate software package.

The archived FASTA is a recovered UniProt reference, not the original search database. The 57 injections were not treated as 57 independent patients. The grouping and reference limitations are described in the manuscript and Supporting Information.

## Using the records

For the exact analysed source, install the pinned commit with `remotes::install_github("enverozturk/atelente", ref = "b20990f7243c87fa468220c7ad27a5af2ce0559a", dependencies = NA)` after installing the dependencies in the [README](../README.md). Repository permission is currently required. The [usage guide](USAGE.md#use-the-installed-package-from-r) describes loading a configuration and calling `atelente::pm_run` with an explicit `config_base`. Use the matching archived input base, adjust only filesystem locations, retain the analytical settings, and select an empty output directory. Read each archive's recorded configuration and execution notes before using its wrapper.

Archive `MANIFEST_SHA256.csv` files identify the included files. Configurations and run manifests also retain input/reference hashes, the analysed commit and recorded environment information. Public sharing copies must document path-only substitutions and their changed file hashes while preserving unchanged input hashes and all scientific values. Original private execution records are retained separately.

The Atelente code is GPL-3; this does not relicense third-party data. Cite the source studies and accessions above and consult the originating repository's data-use terms. Dataset-specific redistribution terms and the cleaned public archive links remain to be checked before depositing the sharing copies.
