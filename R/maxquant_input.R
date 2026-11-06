# MaxQuant peptides.txt sequence summaries -> explicit proteinGroups.txt entities.
# Mod. peptide IDs are source references, not a modification annotation per quantity.
pm_maxquant_flag <- function(data, column) {
  if (!column %in% names(data)) return(rep(NA, nrow(data)))
  value <- trimws(as.character(data[[column]]))
  blank <- is.na(value) | !nzchar(value)
  if (any(!blank & value != "+"))
    stop("MaxQuant ", column, " must contain only '+' or blank values.")
  !blank
}

pm_import_maxquant_tables <- function(peptides, proteins, samples, protein_quantity = "Intensity",
                                      exclude_only_identified_by_site = TRUE) {
  protein_quantity <- pm_input_scalar(protein_quantity, "MaxQuant protein_quantity",
    c("Intensity", "LFQ intensity", "iBAQ"))
  if (!is.logical(exclude_only_identified_by_site) || length(exclude_only_identified_by_site) != 1L ||
      is.na(exclude_only_identified_by_site)) stop("exclude_only_identified_by_site must be TRUE or FALSE.")
  spec <- list(observation_level = "peptide_observation", link_separator = ";", trim_link_whitespace = FALSE,
    input_format = "maxquant", source_grain = "sequence_summary",
    exclude_only_identified_by_site = exclude_only_identified_by_site,
    peptide = list(id = "id", sequence = "Sequence", protein_links = "Protein group IDs",
      scale = "linear", quantity_type = "MaxQuant sequence-summary Intensity"),
    protein = list(id = "id", accessions = "Protein IDs", entity_type = "protein_group",
      scale = "linear", quantity_type = paste("MaxQuant protein-group", protein_quantity)))
  data <- pm_import_quant_tables(peptides, proteins, spec, samples)
  peptide_prefix <- "Intensity "
  protein_prefix <- paste0(protein_quantity, " ")
  peptide_columns <- data$samples$peptide_column
  protein_columns <- data$samples$protein_column
  peptide_experiment <- substring(peptide_columns, nchar(peptide_prefix) + 1L)
  protein_experiment <- substring(protein_columns, nchar(protein_prefix) + 1L)
  if (any(!startsWith(peptide_columns, peptide_prefix) | !nzchar(trimws(peptide_experiment))))
    stop("Select per-experiment peptide 'Intensity <experiment>' columns; total Intensity is not a sample.")
  if (any(!startsWith(protein_columns, protein_prefix) | !nzchar(trimws(protein_experiment))))
    stop("Select per-experiment protein '", protein_prefix, "<experiment>' columns matching protein_quantity.")
  if (any(peptide_experiment != protein_experiment))
    stop("Peptide and protein quantity columns must have exactly matching MaxQuant experiment suffixes.")
  data$samples$source_experiment <- peptide_experiment
  data$observations$source_grain <- rep("sequence_summary", nrow(data$observations))
  data$observations$contaminant_flag <- pm_maxquant_flag(peptides, "Potential contaminant")
  data$observations$decoy_flag <- pm_maxquant_flag(peptides, "Reverse")
  data$entities$contaminant_flag <- pm_maxquant_flag(proteins, "Potential contaminant")
  data$entities$decoy_flag <- pm_maxquant_flag(proteins, "Reverse")
  data$entities$only_identified_by_site_flag <- pm_maxquant_flag(proteins, "Only identified by site")
  append_reason <- function(table, selected, reason) {
    old <- table$exclusion_reason[selected]
    table$exclusion_reason[selected] <- ifelse(nzchar(old), paste(old, reason, sep = ";"), reason)
    table$excluded <- nzchar(table$exclusion_reason)
    table
  }
  for (kind in c("contaminant", "decoy")) {
    flag <- paste0(kind, "_flag")
    data$observations <- append_reason(data$observations, data$observations[[flag]] %in% TRUE,
      paste0("source_", kind, "_flag"))
    data$entities <- append_reason(data$entities, data$entities[[flag]] %in% TRUE,
      paste0("source_", kind, "_flag"))
    excluded_ids <- data$entities$entity_id[data$entities[[flag]] %in% TRUE]
    linked_observations <- unique(data$links$observation_id[data$links$entity_id %in% excluded_ids])
    data$observations <- append_reason(data$observations,
      data$observations$observation_id %in% linked_observations, paste0("linked_", kind, "_entity"))
  }
  if (exclude_only_identified_by_site)
    data$entities <- append_reason(data$entities, data$entities$only_identified_by_site_flag %in% TRUE,
      "only_identified_by_site")
  data
}
