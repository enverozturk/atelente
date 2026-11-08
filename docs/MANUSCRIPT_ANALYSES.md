# Manuscript analyses

The manuscript analyses used Atelente commit [`b20990f7243c87fa468220c7ad27a5af2ce0559a`](https://github.com/enverozturk/atelente/tree/b20990f7243c87fa468220c7ad27a5af2ce0559a). Version 0.1.0 is a release candidate: its runtime files in `R/` and `inst/` are unchanged from that commit. The version update does not change or rerun the reported analyses. Cite the analysed commit when referring to the manuscript results.

**Archive status:** the Drive links below identify the current working archives and may require permission. They are not permanent public deposits. Public archive links, the release date and a persistent release identifier remain to be completed before publication. Archived wrappers preserve the original execution record and can require local path adaptation; this document does not claim they execute unchanged after downloading.

## Data and recorded analyses

All three native analyses used `cleavage.analysis_scope = "specificity_v3"`. Full column mappings, reference settings and sample designs are retained in the configurations below. Paths are relative to the root of the named ZIP, not this repository.

| Study and source | Design used in the manuscript | Working archive |
| --- | --- | --- |
| Schlage et al., 2014, *Time-resolved Analysis of the Matrix Metalloproteinase 10 Substrate Degradome*. [DOI: 10.1074/mcp.M113.035139](https://doi.org/10.1074/mcp.M113.035139); [PXD000503](https://www.ebi.ac.uk/pride/archive/projects/PXD000503). | The publication's **GluC** Supplementary Table S1 time course; pre-digestion iTRAQ, lysine-blocked digestion; one observation per time point, descriptive 16 h versus 1 h comparison, testing disabled. The source article title does not identify the treatment of this selected experiment. | [Schlage sharing archive v3](https://drive.google.com/file/d/1kZCmy0KBIpXky5wAcT5n2vCQVJa8GShq/view): `Schlage_GluC_analyses_and_results_b20990f_PUBLIC_CANDIDATE_v3.zip` |
| Cosenza-Contreras et al., 2024, *TermineR: Extracting information on endogenous proteolytic processing from shotgun proteomics data*. [DOI: 10.1002/pmic.202300491](https://doi.org/10.1002/pmic.202300491); MassIVE accession **MSV000094661**. | Mouse kidney, six WT and five KO samples; pre-digestion TMTpro, lysine-blocked digestion; KO versus WT. The prepared input retains complete observed PSM quantity vectors, with no additional normalization or imputation. | [PKD sharing archive v3](https://drive.google.com/file/d/1MuUpRBYyPT4Qk8qLO5UM1X8n7GQTDGAZ/view): `PKD_analyses_and_results_b20990f_PUBLIC_CANDIDATE_v3_core.zip`; also download [required references](https://drive.google.com/file/d/1GGlYYnnFx6JPs_CogjOV6Tqk-6cidJJ8/view): `PKD_analyses_and_results_b20990f_PUBLIC_CANDIDATE_v3_references.zip` |
| Iacobucci et al., 2026, *Proteomic and degradomic signatures of extracellular matrix remodeling in calcific aortic valve stenosis*. [DOI: 10.1186/s10020-026-01499-0](https://doi.org/10.1186/s10020-026-01499-0); [PXD073246](https://www.ebi.ac.uk/pride/archive/projects/PXD073246). | Non-enriched semi-tryptic **DDA** peptide table; nine Control and ten SEVERE QC-supported triplet means, not named donors. Finite positive quantities were averaged before native ratio testing. Inference is conditional on the F-column grouping; columns were matched by exact header. | [Aortic sharing archive v3](https://drive.google.com/file/d/1GCAwPSMsdnd10Sh3X95UO-2jMydb3Cu9/view): `Aortic_analyses_and_results_b20990f_PUBLIC_CANDIDATE_v3_core.zip`; also download [required references](https://drive.google.com/file/d/1afAz8u1i0PdxvGDCwF7VmuNItIm8bb1z/view): `Aortic_analyses_and_results_b20990f_PUBLIC_CANDIDATE_v3_references.zip` |

### Schlage GluC

- Recorded run configuration: `Schlage_GluC/config_used.json`; existing relative-path configuration: `Schlage_GluC/config_portable.json`.
- Inputs: `Schlage_GluC/inputs/Schlage_peptides_from_audit.csv` and `Schlage_GluC/inputs/Schlage_reference.fasta`.
- Original run wrapper: `Schlage_GluC/provenance/run_schlage_frozen_b20990f.R`.
- Native outputs and run manifest: `Schlage_GluC/results/`; reference snapshots and their provenance: `Schlage_GluC/references/`.

The recovered peptide table is a re-export of archived observations, not a claim of byte identity with the publication's original input file. The configuration and validation records document that distinction.

### PKD and the published-effect/reference comparisons

- Recorded native configuration: `PKD/config_used.json`; relative-path configuration: `PKD/config_portable.json`; input base: `PKD/provenance/input_snapshot/`.
- Exact annotation snapshots and provider provenance: `PKD/references/`.
- Original wrapper: `PKD/provenance/run_pkd_frozen_b20990f.R`; native outputs: `PKD/results/`; execution notes: `PKD/RUN_README.txt`.
- Published S2 comparison: `PKD_comparisons/reproducibility/compare_s2_effects.R`.
- Controlled reference comparison: `PKD_comparisons/reproducibility/compare_pkd_references.R`; model and cohort record: `PKD_comparisons/controlled/analysis_manifest.json`.

The first comparison retains the reported S2 effects. The second compares an observed-peptidoform sum (proxy) with a spanning reference on a common feature/sample set and model; it is not an exact reproduction of TermineR's normalization. Original test families and q values are retained in the archived outputs.

### Aortic valve and technical-replicate preparation

- Recorded configuration and native wrapper: `Aortic/config_used.json`, `Aortic/run_analysis.R`; relative-path configuration: `Aortic/config_portable.json`.
- Exact annotation snapshots and provider provenance: `Aortic/references/`.
- Inputs: `Aortic/inputs/subject_means.csv`, `subject_design.xlsx`, `technical_measurement_counts.csv`, `triplet_mapping.csv` and `recovered_reference_20261104.fasta` (all under `Aortic/inputs/`).
- Native outputs: `Aortic/results/`; extracted audit: `Aortic/audit/`; numerical checks: `Aortic/validation/`.
- Pre-package preparation scripts are included in `Aortic/preprocessing/`: `aort_triplet_qc.R`, `prepare_aort_qc_positive.R` and `prepare_aort_inference.R`. They preserve the historical workflow and may require local path adaptation.

The archived FASTA is a recovered UniProt reference, not the original search database. The 57 injections were not treated as 57 independent patients. The grouping and reference limitations are described in the manuscript and Supporting Information.

## Using the records

The v3 sharing archives include the exact MEROPS and UniProt snapshots recorded by the native audits. Snapshot SHA-256 values match those audits. Scientific tables and input quantities are unchanged; metadata-only substitutions are logged and current checksums are regenerated; no analysis was rerun to complete the archives. The portable configurations retain analytical settings and use archive-relative inputs and reference files. Set `config_base` to the extracted dataset directory (`Schlage_GluC`, `PKD` or `Aortic`) and use a fresh output directory. This packaging check is not a claim of a complete execution rerun.

PKD and Aortic are each delivered as two ZIPs because of a file-upload size limit: core plus required references. Extract both into the same parent directory so that the reference files appear inside the dataset directory. The core ZIP alone is incomplete for annotation-enabled execution. Each ZIP has its own manifest.

Additional historical scripts and Figure 2 artifact provenance are in the [additional-script sharing archive v3](https://drive.google.com/file/d/1wK3BXcHIFvD1LPnak5xy47pDoWrio4GJ/view). The delivered Figure 2 PNG/PDF came from the recorded axis-label correction of the existing figure; the corrected R plotting source is retained separately. No quantitative figure data were changed. Original data-preparation and PowerPoint helper scripts are included; their availability and hashes are recorded alongside the figure provenance.


For the exact analysed source, install the pinned commit with `remotes::install_github("enverozturk/atelente", ref = "b20990f7243c87fa468220c7ad27a5af2ce0559a", dependencies = NA)` after installing the dependencies in the [README](../README.md). Repository permission is currently required. The [usage guide](USAGE.md#use-the-installed-package-from-r) describes loading a configuration and calling `atelente::pm_run` with an explicit `config_base`. Use the matching archived input base, adjust only filesystem locations, retain the analytical settings, and select an empty output directory. Read each archive's recorded configuration and execution notes before using its wrapper.

Archive `MANIFEST_SHA256.csv` files identify the included files. Configurations and run manifests also retain input/reference hashes, the analysed commit and recorded environment information. Sharing copies document path-only substitutions and their changed file hashes while preserving unchanged input hashes and all scientific values. Revision records distinguish the original archives, earlier sharing copies and current files; historical manifests identify the versions named in those records. Original private execution records are retained separately.

The Atelente code is GPL-3; this does not relicense third-party data. Cite the source studies and accessions above and consult the originating repository's data-use terms. Data-use terms apply per source file, not as a blanket license for the combined archive. The article, publisher supplementary tables and repository datasets may have different terms. Source-specific redistribution notices and permanent public links remain to be finalized before public deposit.
