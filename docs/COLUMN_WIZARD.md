# Local Shiny workflow

Install and open Atelente as shown in the [README](../README.md). The four steps are Files, Columns, Samples and Review & run.

Choose the peptide table and its matching FASTA. Select the sequence, master/representative accession set and quantity columns explicitly. Optional modification and contaminant fields remain optional. Quantities may be linear or log2, as declared in the setup.

Give each included quantity column a sample name and condition. A comparison is optional. For a comparison, choose treatment and control. A positive ratio log2FC means a higher neo/spanning ratio in treatment. Optional tests assume independent samples. Technical injections must not be entered as separate independent samples.

Confirm the digestion settings. The default is trypsin. Select lysine blocking only when it applies to the preparation, and declare any applicable pre-digestion N-terminal label. These choices affect candidate classification.

Reference annotations are enabled by default. Their first use requires downloads unless suitable cached references or explicit snapshots are available. Review the destination and run the analysis. The default result folder is `atelente_runs/<run>/results/` under the R working directory.

See [Usage](USAGE.md) for input requirements and downloads, [Analysis rules](ANALYSIS_RULES.md) for calculations, and [Reference annotations](REFERENCE_ANNOTATIONS.md) for cache and offline operation.
