# Attach quantities or refresh an existing report

The main Shiny/`pm_run` workflow already generates its report and export bundle. The commands here support existing reports and the separate peptide/protein workflow. They do not replace the main neo/spanning pair model.

## Attach peptide/protein quantities

Use the original mapping report and a normalization export from the same complete, ordered source peptide table. Optional statistics must come from that exact normalization export. Configure `mapping_dir`, `ratios_dir`, optional `statistics_dir` and `output_dir`, following the [report example](../inst/extdata/config/report.example.json).

```sh
Rscript refresh_report.R --config my_report.json
```

Use a new output directory. Omit `statistics_dir` for quantities without tests. The original mapping report must not already contain an attached quantity layer. Keep the source HTML and adjacent audit artifacts together. Both flat and `audit/` artifact locations are supported.

The bridge checks the complete source table and order, then row IDs, sequence and modification identity. The known mapping-export `peptide_id`/`source_row` prefix is accepted only after identity checks. Matching sequence alone is insufficient. Consumed manifests and table hashes are verified. No mapping, ratio or test is refitted by this attachment step.

## Repackage a report

```sh
Rscript tools/refresh_output_bundle.R INPUT_REPORT NEW_DIRECTORY
```

This command regenerates the compact output bundle, including for a report with attached quantities and optional statistics. Existing artifacts are retained and linked under the output audit folder. It neither attaches a second quantity layer nor recalculates tests.

Output is staged in a temporary sibling directory and published after successful completion. A failed attempt cleans its temporary output and leaves the source bundle unchanged. Keep the resulting folder together for report downloads.

Converting a large historical HTML can require substantially more memory than a fresh analysis. If conversion exceeds available memory, rerun from the original inputs with the current version. A completed compact report preserves the original effects and tests while display filters remain adjustable.
