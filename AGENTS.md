# Working on Atelente

Keep the application and documentation in English. Read [Usage](docs/USAGE.md), [Analysis rules](docs/ANALYSIS_RULES.md) and [Reference annotations](docs/REFERENCE_ANNOTATIONS.md) before changing their contracts.

- Preserve original observations, modifications, literal identifiers, source rows and missing-value states. Alternative protein mappings are not additional independent measurements.
- Separate presentation changes from analysis changes. Preserve ratio effects, original test families and representative selection when editing reports.
- Keep production calculations in R. Browser display filters must not refit tests or reselect the fixed representative.
- Use synthetic fixtures for portable tests. Do not commit experimental inputs, generated reports, caches, credentials or local paths.
- Run the relevant R and JavaScript regression scripts after changes. For package changes, build from the repository root and run installed/source equivalence with an explicitly supplied library path.
- Keep generated test/build output outside the source tree where practical. Verify the package archive contents as well as repository exclusions.
- Update the current documentation when a supported contract changes. Distinguish implemented behavior from proposals and scientific interpretation.

The workflows under docs/GENERIC_INPUT.md, docs/QUANTIFICATION.md, docs/STATISTICS.md and docs/REPORT_QUANTITIES.md are separate peptide/protein interfaces. Do not silently substitute them for the main neo/spanning workflow.
