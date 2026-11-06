# UniProt processing annotations

UniProt annotations describe reference processing boundaries alongside observed candidate bonds. Coordinates are accepted only after matching the accession and complete supplied reference sequence. Canonical coordinates are not automatically transferred to isoforms.

Feature types, boundary uncertainty and evidence qualifiers are retained. Missing reference coverage, sequence mismatch and an evaluated boundary with no matching feature are distinct states. An annotation does not create a quantitative observation or change a fitted test.

See [Reference annotations](REFERENCE_ANNOTATIONS.md) for REST access, local caching, supplied snapshots, provenance and unsupported-identifier handling.
