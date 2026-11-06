# Per-sample quantities over explicit observation-to-protein-entity links.
# Protein groups retain one quantity; their accession lists are never expanded.
pm_empty_sample_ratios <- function() data.frame(
  observation_id=character(), source_row=integer(), sequence=character(), modifications=character(),
  modifications_available=logical(), entity_id=character(), entity_type=character(), accessions=character(),
  entity_source_row=integer(), candidate_entity_count=integer(), sample_id=character(), condition=character(),
  peptide_column=character(), protein_column=character(), peptide_scale=character(), protein_scale=character(),
  peptide_input_value=numeric(), protein_input_value=numeric(), peptide_original_value=character(), protein_original_value=character(),
  peptide_log2_value=numeric(), protein_log2_value=numeric(), normalized_ratio=numeric(), log2_ratio=numeric(),
  peptide_value_status=character(), protein_value_status=character(), entity_status=character(),
  observation_excluded=logical(), entity_excluded=logical(), observation_exclusion_reason=character(), entity_exclusion_reason=character(),
  usable_pair=logical(), ratio_status=character(), ratio_representation=character(), stringsAsFactors=FALSE)

pm_sample_ratios <- function(data) {
  requirements <- list(
    observations=c("observation_id","source_row","sequence","modifications","modifications_available","excluded","exclusion_reason"),
    entities=c("entity_id","entity_type","accessions","source_row","excluded","exclusion_reason"),
    links=c("observation_id","entity_id"), samples=c("sample_id","condition","peptide_column","protein_column"),
    peptide_measurements=c("observation_id","sample_id","source_column","input_value","original_value","value_status","log2_value"),
    protein_measurements=c("entity_id","sample_id","source_column","input_value","original_value","value_status","log2_value"))
  for(name in names(requirements)) {
    table <- data[[name]]
    if(!is.data.frame(table) || anyDuplicated(names(table)) || !all(requirements[[name]] %in% names(table)))
      stop("Missing or duplicate canonical fields in ",name,".")
  }
  id <- function(x,label,unique=TRUE) {
    if(!is.character(x) || anyNA(x) || any(!nzchar(trimws(x))) || (unique && anyDuplicated(x)))
      stop(label," must contain nonempty literal text identifiers",if(unique) " without duplicates" else "",".")
  }
  id(data$observations$observation_id,"Observation IDs")
  id(data$entities$entity_id,"Protein entity IDs")
  id(data$samples$sample_id,"Sample IDs")
  if(!nrow(data$samples)) stop("Include at least one sample for per-sample ratios.")
  for(column in c("peptide_column","protein_column")) id(data$samples[[column]],paste("Selected",column))
  for(table in c("observations","entities")) {
    if(!is.logical(data[[table]]$excluded) || anyNA(data[[table]]$excluded)) stop(table," excluded flags must be explicit logical values.")
  }
  if(anyNA(data$entities$entity_type) || any(!data$entities$entity_type %in% c("protein","protein_group")))
    stop("Protein entity types must be protein or protein_group.")
  id(data$links$observation_id,"Link observation IDs",FALSE)
  id(data$links$entity_id,"Link protein entity IDs",FALSE)
  if(anyDuplicated(data$links[,c("observation_id","entity_id"),drop=FALSE])) stop("Duplicate observation-to-entity links are not independent measurements.")
  if(any(!data$links$observation_id %in% data$observations$observation_id)) stop("Links refer to unknown observations.")
  scales <- c(peptide=data$spec$peptide$scale,protein=data$spec$protein$scale)
  if(length(scales)!=2L || !identical(names(scales),c("peptide","protein")) || anyNA(scales) || any(!scales %in% c("linear","log2")))
    stop("Declare linear or log2 scale for both measurement layers.")
  for(layer in c("peptide","protein")) {
    table <- data[[paste0(layer,"_measurements")]]
    key <- if(layer=="peptide") "observation_id" else "entity_id"
    ids <- if(layer=="peptide") data$observations$observation_id else data$entities$entity_id
    id(table[[key]],paste(layer,"measurement row IDs"),FALSE)
    id(table$sample_id,paste(layer,"measurement sample IDs"),FALSE)
    if(any(!table[[key]] %in% ids) || any(!table$sample_id %in% data$samples$sample_id))
      stop(layer," measurements refer to unknown rows or samples.")
    if(anyDuplicated(table[,c(key,"sample_id"),drop=FALSE])) stop("Duplicate ",layer," measurements for a row/sample.")
    expected <- data$samples[[paste0(layer,"_column")]][match(table$sample_id,data$samples$sample_id)]
    if(!is.character(table$source_column) || anyNA(table$source_column) || any(table$source_column!=expected))
      stop(layer," measurement columns do not match the declared sample mapping.")
    if(!is.numeric(table$input_value) || !is.numeric(table$log2_value)) stop(layer," canonical values must be numeric.")
    if(anyNA(table$value_status) || any(!table$value_status %in% c("observed","missing","zero","negative","invalid_numeric_value")))
      stop("Unknown ",layer," measurement value status.")
    observed <- table$value_status=="observed"
    if(any(observed & (!is.finite(table$input_value) | !is.finite(table$log2_value))))
      stop("Observed ",layer," measurements need finite input and log2 values.")
    if(scales[[layer]]=="linear" && any(observed & table$input_value<=0))
      stop("Observed linear quantities must be positive.")
    expected_log <- if(scales[[layer]]=="linear") log2(table$input_value[observed]) else table$input_value[observed]
    actual_log <- table$log2_value[observed]
    if(any(abs(actual_log-expected_log)>1e-12*pmax(1,abs(expected_log))))
      stop(layer," log2 values disagree with the declared input scale.")
  }
  settings <- list(operation="Per-sample peptide quantity divided by linked protein-entity quantity",
    log_operation="peptide_log2_value - protein_log2_value",scales=as.list(scales),
    quantity_types=list(peptide=data$spec$peptide$quantity_type,protein=data$spec$protein$quantity_type),
    links="Each explicitly linked entity is a candidate denominator; group accessions are not independent entities.",
    observations="Candidate ratio rows are not additional independent peptide measurements.",
    missing="No imputation or pseudocount; missing, zero, negative linear and invalid values stay distinct.",
    range="Finite log2 ratios remain usable when their linear representation overflows or underflows; unavailable linear values remain blank.",
    scope="No condition comparison, replicate-type inference, p-values, site aggregation or spanning-policy decision.")
  summarize <- function(ratios) {
    states <- data.frame(ratio_status=character(),rows=integer(),stringsAsFactors=FALSE)
    if(nrow(ratios)) {
      states <- as.data.frame(table(ratios$ratio_status),stringsAsFactors=FALSE)
      names(states) <- c("ratio_status","rows")
    }
    list(observations=nrow(data$observations),protein_entities=nrow(data$entities),candidate_links=nrow(data$links),
      peptide_measurements=nrow(data$peptide_measurements),protein_measurements=nrow(data$protein_measurements),
      ratio_rows=nrow(ratios),usable_pairs=sum(ratios$usable_pair),linear_ratios=sum(is.finite(ratios$normalized_ratio)),status_counts=states)
  }
  if(!nrow(data$observations)) {
    ratios <- pm_empty_sample_ratios()
    return(list(ratios=ratios,summary=summarize(ratios),settings=settings))
  }
  observations <- data$observations;entities <- data$entities;samples <- data$samples
  # Integer lookup indices, rather than user-key concatenation, prevent collisions.
  link_index <- split(seq_len(nrow(data$links)),match(data$links$observation_id,observations$observation_id))
  counts <- vapply(seq_len(nrow(observations)),function(i) length(link_index[[as.character(i)]]),integer(1))
  candidate_observation <- rep(seq_len(nrow(observations)),pmax(1L,counts))
  candidate_entity <- unlist(lapply(seq_len(nrow(observations)),function(i) {
    rows <- link_index[[as.character(i)]]
    if(length(rows)) data$links$entity_id[rows] else NA_character_
  }),use.names=FALSE)
  nc <- length(candidate_entity);ns <- nrow(samples)
  ci <- rep(seq_len(nc),each=ns);si <- rep(seq_len(ns),nc)
  oi <- candidate_observation[ci];ei <- match(candidate_entity[ci],entities$entity_id)
  ratios <- pm_empty_sample_ratios()[rep(NA_integer_,length(ci)),,drop=FALSE]
  ratios$observation_id <- observations$observation_id[oi]
  for(column in c("source_row","sequence","modifications","modifications_available")) ratios[[column]] <- observations[[column]][oi]
  ratios$entity_id <- candidate_entity[ci];ratios$candidate_entity_count <- counts[oi]
  ratios$entity_type <- entities$entity_type[ei];ratios$accessions <- entities$accessions[ei]
  ratios$entity_source_row <- entities$source_row[ei]
  ratios$sample_id <- samples$sample_id[si];ratios$condition <- samples$condition[si]
  ratios$peptide_column <- samples$peptide_column[si];ratios$protein_column <- samples$protein_column[si]
  ratios$peptide_scale <- rep(scales[["peptide"]],length(ci));ratios$protein_scale <- rep(scales[["protein"]],length(ci))
  ratios$observation_excluded <- observations$excluded[oi]
  ratios$entity_excluded <- entities$excluded[ei]
  ratios$observation_exclusion_reason <- observations$exclusion_reason[oi]
  ratios$entity_exclusion_reason <- entities$exclusion_reason[ei]
  ratios$entity_status <- ifelse(is.na(ratios$entity_id),"missing_protein_link",ifelse(is.na(ei),"protein_entity_not_found","matched"))
  for(layer in c("peptide","protein")) {
    measurement <- data[[paste0(layer,"_measurements")]]
    row_indices <- if(layer=="peptide") match(measurement$observation_id,observations$observation_id) else match(measurement$entity_id,entities$entity_id)
    sample_indices <- match(measurement$sample_id,samples$sample_id)
    target <- if(layer=="peptide") oi else ei
    index <- match(paste(target,si,sep=":"),paste(row_indices,sample_indices,sep=":"))
    ratios[[paste0(layer,"_input_value")]] <- measurement$input_value[index]
    ratios[[paste0(layer,"_original_value")]] <- measurement$original_value[index]
    ratios[[paste0(layer,"_log2_value")]] <- measurement$log2_value[index]
    ratios[[paste0(layer,"_value_status")]] <- measurement$value_status[index]
    ratios[[paste0(layer,"_value_status")]][is.na(index)] <- "measurement_not_exported"
  }
  ratios$ratio_status <- ifelse(ratios$observation_excluded,"excluded_observation",
    ifelse(ratios$entity_status!="matched",ratios$entity_status,
    ifelse(ratios$entity_excluded %in% TRUE,"excluded_protein_entity",
    ifelse(ratios$peptide_value_status!="observed",paste0("peptide_",ratios$peptide_value_status),
    ifelse(ratios$protein_value_status!="observed",paste0("protein_",ratios$protein_value_status),"normalized")))))
  ratios$usable_pair <- ratios$ratio_status=="normalized"
  ratios$normalized_ratio <- rep(NA_real_,nrow(ratios));ratios$log2_ratio <- rep(NA_real_,nrow(ratios))
  ratios$ratio_representation <- rep("unavailable",nrow(ratios))
  good <- which(ratios$usable_pair)
  delta <- ratios$peptide_log2_value[good]-ratios$protein_log2_value[good]
  bad_range <- !is.finite(delta)
  ratios$ratio_status[good[bad_range]] <- "numeric_range_exceeded"
  ratios$usable_pair[good[bad_range]] <- FALSE
  good <- good[!bad_range];delta <- delta[!bad_range]
  ratios$log2_ratio[good] <- delta
  linear <- 2^delta
  overflow <- !is.finite(linear);underflow <- !overflow & linear==0
  ratios$ratio_representation[good] <- ifelse(overflow,"log2_only_overflow",ifelse(underflow,"log2_only_underflow","linear_and_log2"))
  linear[overflow | underflow] <- NA_real_
  ratios$normalized_ratio[good] <- linear
  rownames(ratios) <- NULL
  list(ratios=ratios,summary=summarize(ratios),settings=settings)
}
