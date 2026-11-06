# Reference annotations

Atelente uses MEROPS and UniProt as reference annotations alongside the supplied FASTA. Neither database is bundled. Acquisition and caching are implemented in R.

## Download, cache or opt out

Annotations are enabled by default in `atelente::pm_run()` and Shiny. MEROPS is requested for candidate-bond annotation when a cleavage layer exists. UniProt is requested for accessions in retained mappings. The default cache is `tools::R_user_dir("atelente", "cache")`.

The first required request downloads and validates its reference files. Subsequent runs reuse verified entries without contacting a server. A cached UniProt request can cover a subset of its accession set. New accessions can require another download. Records not returned by a request remain unavailable, rather than becoming negative matches.

To run without annotations, deselect them in Shiny or set `annotations = FALSE` in an R configuration (`"annotations": false` in JSON). This is recorded in the output. Combining this opt-out with explicit reference files is rejected. Required network, response-validation or cache-integrity failures stop the run rather than silently disabling annotation.

UniProt requests send accession identifiers, not peptide tables, quantities or sample names. Unsupported identifier formats remain mapped and are reported as not queried. They do not cause an automatic substitution with another accession.

## Reuse an explicit snapshot

In Shiny, enable **Use my own reference files** and provide the snapshot and provenance paths for each supplied provider. In a run configuration, use:

```json
{
  "merops": {
    "snapshot": "reference/merops_snapshot.tsv",
    "provenance": "reference/merops_provenance.json"
  },
  "uniprot": {
    "snapshot": "reference/uniprot_records.json.gz",
    "provenance": "reference/uniprot_provenance.json"
  }
}
```

These are illustrative paths to validated Atelente-format snapshots and their matching provenance, not arbitrary renamed database files. The snippet is added to a complete run configuration. Relative paths resolve against that configuration's directory. A supplied provider overrides automatic retrieval for that provider only.

Snapshots are not silently refreshed or mixed across releases. For a newer source, use a different cache directory (`cache_dir`) or a new validated snapshot. Hash checks verify consistency with provenance, not the authenticity of an independently supplied provenance document.

## Meaning of a match

**MEROPS:** same-substrate records are found using exact accession and P1 coordinate. The P4–P4′ context is separately classified as agreement, mismatch or unknown. An identical complete eight-residue context elsewhere is retained as motif compatibility, not as a record at this substrate's bond. Protease identifiers, names and source evidence are preserved. MEROPS matching does not verify a full substrate sequence.

**UniProt:** feature coordinates are applied only when the exact primary accession and full reference sequence match the supplied FASTA. Canonical features are not automatically inherited by an isoform, and secondary accessions do not trigger coordinate transfer. Exact signal/transit, propeptide, chain/peptide and supported initiator-Met boundaries can provide processing context. Uncertain feature boundaries and evidence qualifiers remain distinguishable.

Neither source changes the supplied FASTA, creates new candidate sites or refits quantitative tests. A missing record, an unavailable lookup and a sequence mismatch are different states. A reference match is contextual evidence, not a determination of the active protease in the experiment.

## Provenance

The report and manifest record access mode, supplied release/date when available, retrieval time and snapshot hashes. Retained provider provenance includes request/source details. The official standalone MEROPS substrate file may not declare a release number, so it is recorded as unknown rather than inferred from a website label. UniProt responses within one snapshot must agree on release metadata.

For reproducibility, retain the run configuration, `run_manifest.json`, reference snapshots/provenance and R session information locally. Reference databases do not need to be uploaded to the software repository.
