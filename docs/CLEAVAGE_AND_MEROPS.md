# Candidate bonds and MEROPS

The current [analysis rules](ANALYSIS_RULES.md) define candidate bonds, terminal roles, spanning evidence and representative selection. For a bond `k | k+1`, a neo-N form begins at `k+1`, a neo-C form ends at `k`, and a spanning form covers both residues. All coordinates refer to the supplied FASTA.

MEROPS annotations retain the reference substrate, coordinate, context, protease identifier and available evidence. Same-substrate/coordinate matching is distinguished from context-only matching. Context disagreement and unavailable context remain explicit states. A reference match supplies context for a candidate, rather than establishing which enzyme acted in the sample.

Acquisition, provenance and the exact matching requirements are described in [Reference annotations](REFERENCE_ANNOTATIONS.md). References are retrieved and cached in R. No database or external Python acquisition tool is required in the package.
