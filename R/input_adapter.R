# Explicit, software-neutral wide-table import. No abundance estimator or rollup.
pm_input_scalar <- function(value, label, choices = NULL, allow_null = FALSE) {
  if (is.null(value) && allow_null) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value)))
    stop(label, " must be one nonempty string.")
  if (!is.null(choices) && !value %in% choices) stop(label, " must be one of: ", paste(choices, collapse = ", "), ".")
  value
}

pm_input_table <- function(data, label) {
  if (!is.data.frame(data)) stop(label, " must be a data frame.")
  if (anyNA(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data)))
    stop(label, " headers must be nonempty and unique.")
  if (any(vapply(data, is.list, logical(1)))) stop(label, " cannot contain list columns.")
  invisible(data)
}

pm_input_column <- function(data, name, label, required = TRUE) {
  if (is.null(name) && !required) return(NULL)
  name <- pm_input_scalar(name, label)
  if (!name %in% names(data)) stop("Unknown ", label, ": ", name, ".")
  name
}

pm_input_ids <- function(x, label) {
  x <- as.character(x)
  if (anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x))
    stop(label, " must be nonempty, literal, unique identifiers.")
  x
}

pm_input_flags <- function(x, n, label) {
  if (is.null(x)) return(rep(NA, n))
  if (!exists("pm_flag", mode = "function")) stop("Source R/mapping.R before importing tables.")
  original <- as.character(x)
  blank <- is.na(original) | !nzchar(trimws(original))
  flags <- pm_flag(original)
  if (any(!blank & is.na(flags))) stop(label, " contains an unrecognized nonempty flag; use explicit true/false values.")
  flags
}

pm_input_values <- function(x, scale) {
  scale <- pm_input_scalar(scale, "Quantity scale", c("linear", "log2"))
  original <- as.character(x)
  blank <- is.na(original) | trimws(original) %in% c("", "NA", "N/A", "n/a")
  value <- suppressWarnings(as.numeric(original))
  status <- rep("observed", length(value))
  status[!is.finite(value)] <- "invalid_numeric_value"
  status[blank] <- "missing"
  if (scale == "linear") {
    status[!blank & is.finite(value) & value == 0] <- "zero"
    status[!blank & is.finite(value) & value < 0] <- "negative"
  }
  value[blank | !is.finite(value)] <- NA_real_
  logvalue <- rep(NA_real_, length(value))
  usable <- status == "observed"
  logvalue[usable] <- if (scale == "linear") log2(value[usable]) else value[usable]
  data.frame(input_value = value, original_value = original, value_status = status,
             log2_value = logvalue, stringsAsFactors = FALSE)
}

pm_input_measurements <- function(data, ids, samples, id_name, column_name, scale) {
  empty <- data.frame(id = character(), sample_id = character(), source_column = character(),
    input_value = numeric(), original_value = character(), value_status = character(),
    log2_value = numeric(), stringsAsFactors = FALSE)
  names(empty)[1] <- id_name
  if (!nrow(data) || !nrow(samples)) return(empty)
  rows <- lapply(seq_len(nrow(samples)), function(i) {
    column <- samples[[column_name]][i]
    out <- data.frame(id = ids, sample_id = rep(samples$sample_id[i], length(ids)),
      source_column = rep(column, length(ids)), pm_input_values(data[[column]], scale),
      stringsAsFactors = FALSE)
    names(out)[1] <- id_name
    out
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

pm_import_quant_tables <- function(peptides, proteins, spec, samples) {
  pm_input_table(peptides, "Peptide table"); pm_input_table(proteins, "Protein table")
  pm_input_table(samples, "Sample table")
  if (!is.list(spec) || !is.list(spec$peptide) || !is.list(spec$protein))
    stop("Specify peptide and protein column mappings explicitly.")
  pm_input_scalar(spec$observation_level, "observation_level", "peptide_observation")
  if (!"link_separator" %in% names(spec)) stop("Specify link_separator explicitly; use NULL for one literal entity ID per cell.")
  separator <- pm_input_scalar(spec$link_separator, "link_separator", allow_null = TRUE)
  trim_links <- if (is.null(spec$trim_link_whitespace)) FALSE else spec$trim_link_whitespace
  if (!is.logical(trim_links) || length(trim_links) != 1L || is.na(trim_links)) stop("trim_link_whitespace must be true or false.")
  for (layer in c("peptide", "protein")) {
    pm_input_scalar(spec[[layer]]$scale, paste(layer, "scale"), c("linear", "log2"))
    pm_input_scalar(spec[[layer]]$quantity_type, paste(layer, "quantity_type"))
  }
  entity_type <- pm_input_scalar(spec$protein$entity_type, "protein entity_type", c("protein", "protein_group"))
  pcols <- lapply(c("id", "sequence", "protein_links", "modifications", "contaminant"), function(key)
    pm_input_column(peptides, spec$peptide[[key]], paste("peptide", key), required = key %in% c("sequence", "protein_links")))
  names(pcols) <- c("id", "sequence", "protein_links", "modifications", "contaminant")
  rcols <- lapply(c("id", "accessions", "contaminant"), function(key)
    pm_input_column(proteins, spec$protein[[key]], paste("protein", key), required = key == "id"))
  names(rcols) <- c("id", "accessions", "contaminant")
  required_samples <- c("Sample", "Condition", "Use", "Peptide column", "Protein column")
  if (!setequal(names(samples), required_samples)) stop("Sample table requires exactly: ", paste(required_samples, collapse = ", "), ".")
  sample_data <- samples[, required_samples, drop = FALSE]
  for (name in names(sample_data)) sample_data[[name]] <- as.character(sample_data[[name]])
  use <- tolower(trimws(sample_data$Use))
  if (anyNA(use) || any(!use %in% c("yes", "no"))) stop("Use must be Yes or No for every sample row.")
  active <- sample_data[use == "yes", , drop = FALSE]
  if (!nrow(active)) stop("Select at least one sample with Use = Yes.")
  sample_ids <- pm_input_ids(active$Sample, "Selected Sample values")
  for (name in c("Peptide column", "Protein column")) {
    columns <- pm_input_ids(active[[name]], paste("Selected", name, "values"))
    source <- if (name == "Peptide column") peptides else proteins
    unknown <- setdiff(columns, names(source))
    if (length(unknown)) stop("Unknown ", name, ": ", paste(unknown, collapse = ", "), ".")
  }
  active$Condition[is.na(active$Condition)] <- ""
  canonical_samples <- data.frame(sample_id = sample_ids, condition = active$Condition,
    peptide_column = active[["Peptide column"]], protein_column = active[["Protein column"]], stringsAsFactors = FALSE)
  observation_ids <- if (is.null(pcols$id)) sprintf("obs_%07d", seq_len(nrow(peptides))) else pm_input_ids(peptides[[pcols$id]], "Observation IDs")
  entity_ids <- pm_input_ids(proteins[[rcols$id]], "Protein entity IDs")
  if (!exists("pm_parse_sequence", mode = "function")) stop("Source R/mapping.R before importing tables.")
  original_sequence <- as.character(peptides[[pcols$sequence]])
  parsed <- lapply(original_sequence, pm_parse_sequence)
  sequence <- vapply(parsed, function(x) x$sequence, character(1))
  sequence_status <- vapply(parsed, function(x) x$status, character(1))
  pflags <- pm_input_flags(if (is.null(pcols$contaminant)) NULL else peptides[[pcols$contaminant]], nrow(peptides), "Peptide contaminant column")
  rflags <- pm_input_flags(if (is.null(rcols$contaminant)) NULL else proteins[[rcols$contaminant]], nrow(proteins), "Protein contaminant column")
  reasons <- ifelse(sequence_status == "parsed", "", sequence_status)
  append_reason <- function(reasons, selected, value) {
    reasons[selected] <- ifelse(nzchar(reasons[selected]), paste(reasons[selected], value, sep = ";"), value)
    reasons
  }
  reasons <- append_reason(reasons, pflags %in% TRUE, "source_contaminant_flag")
  observations <- data.frame(observation_id = observation_ids, source_row = seq_len(nrow(peptides)) + 1L,
    original_sequence = original_sequence, sequence = sequence, sequence_status = sequence_status,
    modifications = if (is.null(pcols$modifications)) rep(NA_character_, nrow(peptides)) else as.character(peptides[[pcols$modifications]]),
    modifications_available = rep(!is.null(pcols$modifications), nrow(peptides)), contaminant_flag = pflags,
    excluded = nzchar(reasons), exclusion_reason = reasons, stringsAsFactors = FALSE)
  entities <- data.frame(entity_id = entity_ids, entity_type = rep(entity_type, nrow(proteins)),
    accessions = if (is.null(rcols$accessions)) rep(NA_character_, nrow(proteins)) else as.character(proteins[[rcols$accessions]]),
    source_row = seq_len(nrow(proteins)) + 1L, contaminant_flag = rflags, excluded = rflags %in% TRUE,
    exclusion_reason = ifelse(rflags %in% TRUE, "source_contaminant_flag", ""), stringsAsFactors = FALSE)
  raw_links <- as.character(peptides[[pcols$protein_links]])
  split_links <- lapply(raw_links, function(value) {
    if (is.na(value) || !nzchar(trimws(value))) return(character())
    tokens <- if (is.null(separator)) value else strsplit(value, separator, fixed = TRUE)[[1]]
    if (trim_links) tokens <- trimws(tokens)
    unique(tokens[!is.na(tokens) & nzchar(trimws(tokens))])
  })
  links <- data.frame(observation_id = rep(observation_ids, lengths(split_links)),
    entity_id = as.character(unlist(split_links, use.names = FALSE)), stringsAsFactors = FALSE)
  contaminated_ids <- entity_ids[rflags %in% TRUE]
  contaminated_observations <- unique(links$observation_id[links$entity_id %in% contaminated_ids])
  observations$exclusion_reason <- append_reason(observations$exclusion_reason,
    observation_ids %in% contaminated_observations, "linked_contaminant_entity")
  observations$excluded <- nzchar(observations$exclusion_reason)
  peptide_measurements <- pm_input_measurements(peptides, observation_ids, canonical_samples,
    "observation_id", "peptide_column", spec$peptide$scale)
  protein_measurements <- pm_input_measurements(proteins, entity_ids, canonical_samples,
    "entity_id", "protein_column", spec$protein$scale)
  list(observations = observations, entities = entities, links = links,
    peptide_measurements = peptide_measurements, protein_measurements = protein_measurements,
    samples = canonical_samples, source_samples = samples, source_peptides = peptides,
    source_proteins = proteins, spec = spec)
}
