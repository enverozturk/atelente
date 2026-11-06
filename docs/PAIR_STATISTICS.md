# Neo/spanning pair statistics

The main workflow can test sample log2 neo/spanning ratios using an independent-sample OLS model with limma moderation. The [analysis rules](ANALYSIS_RULES.md) define eligibility, test families, representative selection and BH adjustment.

Ratios require both forms in the same sample. This is distinct from pairing a control sample with a treatment sample. At least three usable independent samples per condition and at least two eligible pairs in the moderation pool are required. Counts are minimum requirements, not a power guarantee.

A directed measurement pair is tested once per contrast and test family. Alternative protein mappings do not add tests. BH adjustment uses that contrast/family's tested pairs. Selecting a report filter or representative view does not refit the model or redefine the family. Representative selection does not use effect size or p-values.

Missing values are not imputed. A form measured only in one condition receives an availability status instead of an invented finite fold change. Condition-dependent missingness and shared spanning denominators can affect interpretation. Differential ratios do not estimate a cleavage percentage or establish an active protease.

The separate peptide/protein workflow has its own [statistics interface](STATISTICS.md). Its tests and minimum counts must not be substituted for the neo/spanning pair tests.
