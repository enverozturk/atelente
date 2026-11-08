\name{atelente}
\alias{atelente}
\alias{pm_run}
\alias{pm_annotate_references}
\alias{pm_read_table}
\alias{pm_read_fasta}
\alias{pm_map}
\alias{pm_build_cleavage}
\alias{pm_annotate_merops}
\alias{pm_annotate_uniprot}
\alias{pm_neo_spanning_ratios}
\alias{pm_neo_spanning_contrasts}
\alias{pm_read_sample_metadata}
\alias{pm_export}
\alias{pm_import_quant_tables}
\alias{pm_import_maxquant_tables}
\alias{pm_sample_ratios}
\alias{pm_test_ratios}
\alias{pm_test_pair_matrix}
\alias{pm_pair_workflow}
\alias{pm_prepare_column_input}
\alias{pm_wizard_app}
\title{Map terminal and spanning peptides to a common candidate bond}
\description{The primary application is shotgun proteomics without terminal peptide enrichment.
Mapping preserves the original peptide observations, protein candidates and reference coordinates.}
\details{
\code{pm_run(config, config_base = getwd(), template = NULL)} runs the complete
mapping workflow from a named configuration list. Required fields are
\code{peptides}, \code{fasta} and \code{output_dir}; paths resolve relative to
\code{config_base}. The installed package uses its own HTML assets unless
\code{template} is supplied. Add a \code{cleavage} profile for candidate bonds.
Reference annotations default to enabled. Initial access downloads official
references in R and subsequent runs reuse verified local cache entries.
\code{annotations = FALSE} in the configuration explicitly disables annotation.
\code{cache_dir} optionally overrides \code{tools::R_user_dir("atelente", "cache")}.
Each optional \code{merops} or \code{uniprot} list with \code{snapshot} and
\code{provenance} paths overrides automatic access for that source.
\code{pm_annotate_references(result, annotations = TRUE, merops = NULL,
uniprot = NULL, cache_dir = NULL)} exposes the same reference stage for users
composing their own workflow. Required download or cache-verification failures
stop the run with an error. Reference retrieval and annotation run in R.

Use \code{pm_read_table} and \code{pm_read_fasta} to read the experiment and reference sequence.
Pass the mapped result, the declared digestion profile and the reference file SHA-256 to
\code{pm_build_cleavage}. The result contains sites, evidence and annotated peptide termini.
\code{pm_annotate_merops} and \code{pm_annotate_uniprot} accept local snapshots and their provenance files.

\code{pm_neo_spanning_ratios} calculates observed ratios for a terminal peptide and a spanning
peptide measured in the same sample. Use \code{pm_read_sample_metadata} and
\code{pm_neo_spanning_contrasts} for descriptive comparisons.
\code{pm_export} writes the source audits, result tables and a local HTML report to an empty directory.
The CLI and Shiny request \code{compact_bundle=TRUE}, producing a three-sheet
\code{atelente_results.xlsx}, optional offline worksheet selection, and a separate
\code{atelente_audit.zip}. Low-level calls keep the legacy layout unless requested.
Its optional \code{representative_selection_k} argument freezes the positive integer
common-sample minimum used to select report representatives; by default this is
the smaller of three and the smaller included condition size. Display coverage
filters do not reselect representatives. Original pair tests and BH families remain unchanged.

The separate quantity workflow uses \code{pm_import_quant_tables} or
\code{pm_import_maxquant_tables} with explicit sample and protein-group links.
\code{pm_sample_ratios} calculates peptide/protein ratios and \code{pm_test_ratios} compares conditions.
\code{pm_test_pair_matrix} accepts a log2 neo/spanning matrix, independent sample metadata,
pair-family membership and a directional contrast. It tests measurement pairs within each declared
family. Its results describe differential ratios at candidate bonds.

\code{pm_wizard_app()} opens a local Shiny application for explicit column,
quantity-scale, sample and reference selection. Shiny is optional for scripted
use. \code{pm_prepare_column_input} preserves the original table while preparing
the selected columns for mapping. Saved configurations can be passed to
\code{pm_run}; individual sample conditions and independent units are explicit.

The opt-in \code{analysis_scope = "specificity_v3"} profile excludes unspecific
placements from candidate evidence and preserves a separate exclusion audit.
\code{pm_pair_workflow} adapts exact mapping links to the independent-sample
pair model. Representative assignment requires a selected source column with
one common accession and unambiguous placements for both forms. Alternative
mapping links do not add tests. Tests are refitted within each family/contrast
when this scope is used; report filters never recalculate BH.
}
\examples{
ref <- tempfile(fileext = ".fasta")
writeLines(c(">P1", "MACDEFGHIK"), ref)
p <- data.frame(Sequence = c("MACD", "EFGHIK", "CDEFGH"),
                "Master Protein Accessions" = "P1", check.names = FALSE)
x <- pm_map(p, pm_read_fasta(ref), sequence_column = "Sequence")
y <- pm_build_cleavage(x, list(enzyme = "trypsin", lysine_blocked = FALSE),
                       digest::digest(file = ref, algo = "sha256"))
y$cleavage$sites[, c("bond_left", "neo_n_count", "neo_c_count", "spanning_count")]
}
