# Atelente 0.1.0 (unreleased)

- Map peptide observations to reference sequences and link neo-terminal and spanning forms at candidate bonds.
- Calculate sample ratios and condition effects, with optional limma tests and BH adjustment.
- Analyse data through a local Shiny interface and export HTML reports, Excel tables and audit records.
- Annotate candidate bonds with MEROPS and UniProt using downloaded or supplied reference snapshots.
- Group source-checkout commands under `scripts/` and consolidate the usage and method guides.
- Report `Spanning ends` and `Spanning terminal form` in Results and All pairs. Selection, ratios and tests are unchanged.
- Correct the analysis guide: reference ends, initiator-Met boundaries and labelled N termini do not satisfy the end-compatibility preference, as implemented.
- Add `scripts/compare_results.R` to confirm that two result workbooks differ only by added columns.
- Add the three manuscript datasets, saved settings and a shared runner; update usage and citation metadata.
