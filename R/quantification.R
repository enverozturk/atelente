# Descriptive peptide effects relative to supplied PD protein abundances.
# Alternative mapped PD masters are separate candidate baselines, not repeats.
pm_quant_flag <- function(x) {
  if (is.logical(x)) return(x)
  x <- trimws(as.character(x))
  out <- rep(NA, length(x))
  out[x %in% c("TRUE", "True", "true", "1", "Yes", "yes", "YES", "DO\u011eRU")] <- TRUE
  out[x %in% c("FALSE", "False", "false", "0", "No", "no", "NO", "YANLI\u015e")] <- FALSE
  out
}

pm_quant_values <- function(x) {
  original <- as.character(x)
  blank <- is.na(original) | trimws(original) %in% c("", "NA", "N/A", "n/a")
  value <- suppressWarnings(as.numeric(original))
  status <- ifelse(blank, "missing", ifelse(!is.finite(value), "invalid_numeric_value",
    ifelse(value == 0, "zero", ifelse(value < 0, "negative", "positive"))))
  value[blank | !is.finite(value)] <- NA_real_
  list(value = value, status = status, original = original)
}

pm_empty_quant_effects <- function() data.frame(comparison = character(), treatment = character(), control = character(),
  peptide_id = character(), source_row = integer(), sequence = character(), modifications = character(),
  accession = character(), candidate_accession_count = integer(), protein_source_row = integer(),
  protein_identified_peptide_count = numeric(), protein_psm_count = numeric(), protein_quantifying_peptide_count = numeric(),
  raw_peptide_log2fc = numeric(), matched_peptide_log2fc = numeric(), protein_log2fc = numeric(), adjusted_log2fc = numeric(),
  n_peptide_treatment = integer(), n_peptide_control = integer(), n_paired_treatment = integer(), n_paired_control = integer(),
  n_selected_treatment = integer(), n_selected_control = integer(), method = character(),
  protein_reference_status = character(), status = character(), stringsAsFactors = FALSE)

pm_empty_quant_pairs <- function() data.frame(peptide_id = character(), accession = character(),
  candidate_accession_count = integer(), protein_source_row = integer(), sample = character(), condition = character(),
  peptide_column = character(), protein_column = character(), peptide_value = numeric(), protein_value = numeric(),
  peptide_original_value = character(), protein_original_value = character(), peptide_value_status = character(),
  protein_value_status = character(), protein_reference_status = character(), usable_pair = logical(), stringsAsFactors = FALSE)

pm_adjust_peptide_effects <- function(result, metadata, method = "mean_log2") {
  if (!identical(method, "mean_log2")) stop("Only the explicit descriptive method 'mean_log2' is supported.")
  if (!exists("pm_validate_sample_metadata", mode = "function")) stop("Source R/sample_metadata.R first.")
  if (is.null(result$source_proteins) || !"Accession" %in% names(result$source_proteins))
    stop("Protein correction requires a supplied protein table with an Accession column.")
  metadata <- pm_validate_sample_metadata(metadata, result, require_protein = TRUE)
  samples <- metadata$samples[metadata$samples$Use == "Yes", , drop = FALSE]
  comparisons <- metadata$comparisons
  peptides <- result$peptides; matches <- result$matches; measurements <- result$measurements
  required <- list(peptides = c("peptide_id", "mapping_status"),
    matches = c("peptide_id", "accession", "is_master"),
    measurements = c("peptide_id", "source_column", "value", "original_value", "value_status"))
  for (name in names(required)) if (is.null(result[[name]]) || !all(required[[name]] %in% names(result[[name]])))
    stop("Missing required quantification fields in ", name, ".")
  if (anyNA(peptides$peptide_id) || any(!nzchar(peptides$peptide_id)) || anyDuplicated(peptides$peptide_id))
    stop("Peptide IDs must be nonempty and unique.")
  if (any(!matches$peptide_id %in% peptides$peptide_id)) stop("Placements refer to unknown peptide IDs.")
  measurements <- measurements[measurements$source_column %in% samples[["Peptide column"]], , drop = FALSE]
  if (any(!measurements$peptide_id %in% peptides$peptide_id)) stop("Measurements refer to unknown peptide IDs.")
  if (anyDuplicated(measurements[, c("peptide_id", "source_column"), drop = FALSE]))
    stop("A peptide/source-column measurement is duplicated; do not count repeated exports as new observations.")
  if (anyNA(measurements$value_status) || any(!measurements$value_status %in% c("observed", "missing", "invalid_numeric_value")))
    stop("Unknown peptide measurement value status.")
  settings <- list(method = method, direction = "Treatment minus Control on the log2 scale",
    formula = "matched_peptide_log2fc - protein_log2fc",
    sample_rule = "Within each condition, adjusted effects use the same positive peptide/protein sample pairs.",
    baseline_rule = "Each distinct mapped PD master accession is a separate exact-accession candidate baseline.",
    multiplicity = "Candidate baseline rows are not independent peptide observations or measurements.",
    raw_effect = "raw_peptide_log2fc uses all positive peptide measurements; matched_peptide_log2fc uses only usable protein pairs.",
    protein_abundance = "Supplied calculated protein abundance; no recalculation, site-specific normalization or purity assumption.",
    support = "Identified peptide and PSM counts do not establish the peptide support used for quantification.",
    missing = "Missing, invalid, zero and negative values remain distinct; no imputation or pseudocount.",
    inference = "Descriptive effects only. No p-values, confidence intervals, or significance thresholds.",
    site_scope = "No site aggregation or decision about applying protein correction to spanning evidence.")
  keep <- peptides$peptide_id %in% matches$peptide_id &
    !is.na(peptides$mapping_status) & peptides$mapping_status %in% c("one_match_in_scope", "multiple_matches_in_scope")
  active <- peptides[keep, , drop = FALSE]
  if (!nrow(active)) return(list(effects = pm_empty_quant_effects(), measurement_pairs = pm_empty_quant_pairs(), settings = settings))

  masters <- unique(matches[pm_quant_flag(matches$is_master) %in% TRUE &
    matches$peptide_id %in% active$peptide_id & !is.na(matches$accession) & nzchar(matches$accession), c("peptide_id", "accession"), drop = FALSE])
  master_index <- split(masters$accession, masters$peptide_id)
  counts <- vapply(active$peptide_id, function(id) length(master_index[[id]]), integer(1))
  candidate_peptide_index <- rep(seq_len(nrow(active)), pmax(counts, 1L))
  candidate_accessions <- unlist(lapply(active$peptide_id, function(id) {
    accession <- master_index[[id]]
    if (length(accession)) as.character(accession) else NA_character_
  }), use.names = FALSE)
  nc <- length(candidate_accessions); ns <- nrow(samples)
  candidate_ids <- active$peptide_id[candidate_peptide_index]
  candidate_counts <- counts[candidate_peptide_index]
  proteins <- result$source_proteins
  protein_index <- split(seq_len(nrow(proteins)), as.character(proteins$Accession))
  protein_rows <- vapply(candidate_accessions, function(accession) {
    if (is.na(accession)) return(NA_integer_)
    idx <- protein_index[[accession]]
    if (length(idx) == 1L) as.integer(idx) else NA_integer_
  }, integer(1))
  reference <- vapply(candidate_accessions, function(accession) {
    if (is.na(accession)) return("no_mapped_pd_master")
    n <- length(protein_index[[accession]])
    if (!n) "protein_accession_not_found" else if (n > 1L) "ambiguous_protein_rows" else "matched"
  }, character(1))
  if ("Contaminant" %in% names(proteins)) {
    flagged <- pm_quant_flag(proteins$Contaminant)[protein_rows] %in% TRUE
    reference[flagged] <- "protein_flagged_contaminant"
  }

  # Construct a candidate/sample audit once, not once per requested comparison.
  ci <- rep(seq_len(nc), each = ns); si <- rep(seq_len(ns), nc)
  pairs <- pm_empty_quant_pairs()[rep(NA_integer_, length(ci)), , drop = FALSE]
  pairs$peptide_id <- candidate_ids[ci]; pairs$accession <- candidate_accessions[ci]
  pairs$candidate_accession_count <- candidate_counts[ci]
  pairs$protein_source_row <- protein_rows[ci] + 1L
  pairs$sample <- samples$Sample[si]; pairs$condition <- samples$Condition[si]
  pairs$peptide_column <- samples[["Peptide column"]][si]; pairs$protein_column <- samples[["Protein column"]][si]
  pairs$protein_reference_status <- reference[ci]
  # Integer source-column indices prevent collisions even for unusual literal labels.
  measure_column_index <- match(measurements$source_column, samples[["Peptide column"]])
  measure_peptide_index <- match(measurements$peptide_id, peptides$peptide_id)
  keys <- paste(measure_peptide_index, measure_column_index, sep = ":")
  target_keys <- paste(match(pairs$peptide_id, peptides$peptide_id), si, sep = ":")
  mi <- match(target_keys, keys)
  pval <- pm_quant_values(measurements$value[mi])
  pairs$peptide_value <- pval$value
  pairs$peptide_original_value <- as.character(measurements$original_value[mi])
  pairs$peptide_value_status <- pval$status
  source_status <- measurements$value_status[mi]
  pairs$peptide_value_status[source_status %in% c("missing", "invalid_numeric_value")] <- source_status[source_status %in% c("missing", "invalid_numeric_value")]
  pairs$peptide_value[!source_status %in% "observed"] <- NA_real_
  pairs$peptide_value_status[is.na(mi)] <- "measurement_not_exported"
  if ("included_in_analysis" %in% names(measurements)) {
    excluded <- !is.na(mi) & !pm_quant_flag(measurements$included_in_analysis[mi]) %in% TRUE
    pairs$peptide_value_status[excluded] <- "excluded_measurement"
  }
  pairs$protein_value <- rep(NA_real_, nrow(pairs))
  pairs$protein_original_value <- rep(NA_character_, nrow(pairs))
  pairs$protein_value_status <- rep("protein_reference_unavailable", nrow(pairs))
  for (s in seq_len(ns)) {
    rows <- which(si == s & !is.na(protein_rows[ci]))
    values <- pm_quant_values(proteins[[samples[["Protein column"]][s]]][protein_rows[ci[rows]]])
    pairs$protein_value[rows] <- values$value
    pairs$protein_original_value[rows] <- values$original
    pairs$protein_value_status[rows] <- values$status
  }
  pairs$usable_pair <- pairs$peptide_value_status == "positive" & pairs$protein_value_status == "positive" &
    pairs$protein_reference_status == "matched"
  rownames(pairs) <- NULL

  # The first dimension is sample, so each candidate column stays aligned.
  peptide_values <- matrix(pairs$peptide_value, nrow = ns)
  protein_values <- matrix(pairs$protein_value, nrow = ns)
  positive_peptide <- matrix(pairs$peptide_value_status == "positive", nrow = ns)
  paired <- matrix(pairs$usable_pair, nrow = ns)
  logmean <- function(values, valid, selected) {
    used <- valid & matrix(selected, nrow = ns, ncol = nc)
    logged <- matrix(0, nrow = ns, ncol = nc)
    logged[used] <- log2(values[used])
    n <- colSums(used)
    value <- colSums(logged) / n
    value[n == 0L] <- NA_real_
    list(value = value, n = as.integer(n))
  }
  support <- function(field) {
    if (!field %in% names(proteins)) return(rep(NA_real_, nc))
    value <- suppressWarnings(as.numeric(as.character(proteins[[field]][protein_rows])))
    value[!is.finite(value) | value < 0 | value != floor(value)] <- NA_real_
    value
  }
  result_rows <- vector("list", nrow(comparisons))
  for (i in seq_len(nrow(comparisons))) {
    treatment <- samples$Condition == comparisons$Treatment[i]
    control <- samples$Condition == comparisons$Control[i]
    raw_t <- logmean(peptide_values, positive_peptide, treatment); raw_c <- logmean(peptide_values, positive_peptide, control)
    pep_t <- logmean(peptide_values, paired, treatment); pep_c <- logmean(peptide_values, paired, control)
    pro_t <- logmean(protein_values, paired, treatment); pro_c <- logmean(protein_values, paired, control)
    effects <- pm_empty_quant_effects()[rep(NA_integer_, nc), , drop = FALSE]
    effects$comparison <- rep(comparisons$Name[i], nc); effects$treatment <- rep(comparisons$Treatment[i], nc)
    effects$control <- rep(comparisons$Control[i], nc); effects$peptide_id <- candidate_ids
    for (field in c("source_row", "sequence", "modifications")) if (field %in% names(active))
      effects[[field]] <- active[[field]][candidate_peptide_index]
    effects$accession <- candidate_accessions; effects$candidate_accession_count <- candidate_counts
    effects$protein_source_row <- protein_rows + 1L
    effects$protein_identified_peptide_count <- support("# Peptides")
    effects$protein_psm_count <- support("# PSMs")
    # None of the inspected source headers establishes actual quantifying support.
    effects$protein_quantifying_peptide_count <- rep(NA_real_, nc)
    effects$raw_peptide_log2fc <- raw_t$value - raw_c$value
    effects$matched_peptide_log2fc <- pep_t$value - pep_c$value
    effects$protein_log2fc <- pro_t$value - pro_c$value
    effects$adjusted_log2fc <- effects$matched_peptide_log2fc - effects$protein_log2fc
    effects$n_peptide_treatment <- raw_t$n; effects$n_peptide_control <- raw_c$n
    effects$n_paired_treatment <- pep_t$n; effects$n_paired_control <- pep_c$n
    effects$n_selected_treatment <- rep(sum(treatment), nc); effects$n_selected_control <- rep(sum(control), nc)
    effects$method <- rep(method, nc); effects$protein_reference_status <- reference
    effects$status <- ifelse(reference != "matched", reference,
      ifelse(pep_t$n == 0 & pep_c$n == 0, "no_usable_pairs_in_either_condition",
      ifelse(pep_t$n == 0, "no_usable_treatment_pairs", ifelse(pep_c$n == 0, "no_usable_control_pairs", "adjusted"))))
    result_rows[[i]] <- effects
  }
  effects <- do.call(rbind, result_rows); rownames(effects) <- NULL
  list(effects = effects, measurement_pairs = pairs, settings = settings)
}
