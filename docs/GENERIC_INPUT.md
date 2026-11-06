# Explicit peptide/protein quantity input

`normalize.R` is a separate workflow for explicitly linked peptide and protein/group quantities. It calculates per-sample peptide/protein ratios without FASTA placement. For the main neo/spanning workflow, use [Usage](USAGE.md).

## Inputs and configuration

Supply peptide and protein tables, a sample table and a configuration. CSV accepts comma or semicolon delimiters. TSV/TXT and XLSX are also supported. Selected numeric text uses decimal points without thousands separators. Malformed selected quantities stop the run with source row/column details.

A minimal explicit specification is:

```json
{
  "peptides": "peptides.csv",
  "proteins": "proteins.csv",
  "samples": "samples.csv",
  "output_dir": "ratios_results",
  "excel": true,
  "spec": {
    "observation_level": "peptide_observation",
    "link_separator": ";",
    "peptide": {
      "id": "Observation",
      "sequence": "Sequence",
      "protein_links": "Protein groups",
      "modifications": "Modification",
      "scale": "linear",
      "quantity_type": "Reported peptide intensity"
    },
    "protein": {
      "id": "Group",
      "accessions": "Members",
      "entity_type": "protein_group",
      "scale": "linear",
      "quantity_type": "Reported protein-group abundance"
    }
  }
}
```

Use exact source headers. Peptide sequence/links, protein ID, both scales and quantity-type descriptions are required. Peptide ID/modifications/contaminant and protein accessions/contaminant mappings are optional. Missing peptide IDs are generated for that import. Original rows and modification variants remain distinct.

Protein IDs must be nonempty and unique. `entity_type` is `protein` or `protein_group`. `link_separator` separates candidate IDs, not group membership. Use `null` for a single literal ID. Whitespace is preserved unless `trim_link_whitespace` is explicitly true. Repeated links within an observation are deduplicated. Unknown links remain diagnostic.

The Samples table has exactly these columns:

| Sample | Condition | Use | Peptide column | Protein column |
| --- | --- | --- | --- | --- |
| Control_1 | Control | Yes | Intensity 1 | Abundance 1 |
| Treatment_1 | Treatment | Yes | Intensity 2 | Abundance 2 |

Sample names and included quantity columns must be unique. `Use` is `Yes` or `No`. Conditions may be blank or identical. A same-condition run needs no invented comparison. Worksheets can be selected with `peptide_sheet`, `protein_sheet` and `sample_sheet`.

## Calculation and output

```sh
Rscript normalize.R --config my_normalize.json
```

Paths resolve relative to the configuration. Use a new output directory. `scale` is `linear` or `log2` for each layer independently. Positive finite linear values and finite log2 values are usable. The calculation is `peptide_log2 - protein_log2`. No imputation or pseudocount is added. Linear ratio underflow/overflow is flagged while a finite log2 ratio remains available.

Optional explicit contaminant flags exclude the corresponding observations. A contaminated linked entity excludes its observation conservatively. Original flags and reasons are retained. Multiple candidate entities yield alternative denominators, not independent measurements.

Exports include `sample_ratios.csv`, original and canonical tables, `normalization_manifest.json`, and optionally `sample_ratios.xlsx`. IDs remain text and measurements remain numeric. The manifest records inputs, quantities, settings and hashes. The command does not produce an HTML report.

The installed R API exposes `pm_import_quant_tables(peptides, proteins, spec, samples)` followed by `pm_sample_ratios(data)`. Canonical tables retain observations, entities, literal links, measurements, samples and the original inputs.

Optional downstream steps are [statistics](STATISTICS.md) and [attachment to a matching mapping report](REPORT_QUANTITIES.md). See [MaxQuant input](MAXQUANT_VALIDATION.md) for its separate explicit preset.
