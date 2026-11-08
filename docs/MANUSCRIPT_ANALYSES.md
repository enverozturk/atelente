# Manuscript data and analyses

The manuscript analyses used commit [b20990f](https://github.com/enverozturk/atelente/tree/b20990f7243c87fa468220c7ad27a5af2ce0559a). The 0.1.0 candidate retains the same computational code.

| Dataset | Source | Analysis |
| --- | --- | --- |
| Schlage | [PXD000503](https://www.ebi.ac.uk/pride/archive/projects/PXD000503); [paper](https://doi.org/10.1074/mcp.M113.035139) | Published GluC time course; descriptive, no statistical testing |
| PKD | MSV000094661; [TermineR paper](https://doi.org/10.1002/pmic.202300491) | Six WT and five KO samples |
| Aortic valve | [PXD073246](https://www.ebi.ac.uk/pride/archive/projects/PXD073246); [paper](https://doi.org/10.1186/s10020-026-01499-0) | Nine Control and ten SEVERE QC-supported triplet means; not named donors |

The [analysis archive index](https://drive.google.com/file/d/1a0HeQokmYq-FII6ryM4Ctrw24Z6FDpXq/view) links the inputs, FASTA files, configurations, reference snapshots, scripts and results. Its README files describe the recorded analyses and file checksums. These are private working links pending permanent public deposition.

For PKD and aortic valve, extract both the core and reference ZIPs into the same parent directory. Use the dataset's `config_portable.json` with `config_base` set to that dataset directory; see [scripted usage](USAGE.md#use-the-installed-package-from-r). Historical wrappers may require local path adaptation.

Third-party data retain their source terms; the package's GPL-3 license does not relicense them.
