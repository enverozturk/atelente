# Reference annotations

Atelente adds MEROPS substrate records and UniProt processing features to candidate bonds mapped on the supplied FASTA.

## Downloads and offline use

Annotations are enabled by default in Shiny and `pm_run()`. MEROPS is requested for candidate bonds; UniProt is requested for accessions in retained mappings. The default cache is `tools::R_user_dir("atelente", "cache")`.

Downloaded files are checked and cached. Subsequent runs reuse matching entries offline. A cached UniProt response can cover a subset of its accession set; additional accessions may require a download. UniProt requests send accession identifiers. Input tables and sample quantities stay local.

To disable annotations, deselect them in Shiny or set `"annotations": false` in the run configuration. Use this setting without explicit reference-file entries. Download, response-validation and cache-integrity failures produce an error with the required next step.

## Supplied reference snapshots

In Shiny, select **Use my own reference files** and provide each snapshot with its provenance file. Equivalent configuration fields are:

```json
{
  "merops": {
    "snapshot": "references/merops_substrates.csv.gz",
    "provenance": "references/merops_provenance.json"
  },
  "uniprot": {
    "snapshot": "references/uniprot_records.json.gz",
    "provenance": "references/uniprot_provenance.json"
  }
}
```

Add these fields to the run configuration using Atelente-format snapshots and their matching provenance files. Paths resolve relative to that configuration. A supplied snapshot replaces automatic retrieval for that provider.

To use a newer reference version, choose a new `cache_dir` or supply a new snapshot. SHA-256 checks compare snapshot contents with the recorded provenance.

## Matching rules

**MEROPS:** same-substrate matches require the exact accession and P1 coordinate. P4–P4′ context is reported as agreement, mismatch or unavailable. A complete matching eight-residue context at another substrate location is recorded as motif compatibility. Results retain protease identifiers, names and source evidence.

**UniProt:** feature coordinates require an exact primary accession and full-sequence match to the supplied FASTA. This includes signal/transit peptides, propeptides, chains/peptides and supported initiator-Met boundaries. Feature uncertainty and evidence qualifiers remain visible. Unsupported identifiers, missing records and sequence mismatches receive separate statuses.

Annotations provide reference context alongside the measured peptide ratios. Candidate generation and quantitative testing use the experimental input and declared analysis settings.

## Recorded provenance

The report and manifest record access mode, retrieval time, available release metadata and snapshot hashes. Provider records retain source and request details. Unspecified release numbers remain unavailable; UniProt responses within a snapshot must have consistent release metadata.

Keep the configuration, reference snapshots, provenance files and run manifest with the analysis results.
