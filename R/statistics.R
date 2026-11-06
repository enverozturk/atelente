# Moderated two-group tests of per-sample log2 peptide/protein ratios.
# Each selected sample declares an independent experimental unit explicitly.
pm_test_ratios <- function(ratios, samples, comparisons, design="independent_units") {
  if(!is.character(design) || length(design)!=1L || is.na(design) || design!="independent_units")
    stop("Only design = independent_units is implemented; paired, repeated or batch designs require a separate model.")
  exact_table <- function(x,columns,label) {
    if(!is.data.frame(x) || anyDuplicated(names(x)) || !setequal(names(x),columns))
      stop(label," requires exactly: ",paste(columns,collapse=", "),".")
  }
  text_ids <- function(x,label,unique=FALSE) {
    if(!is.character(x) || anyNA(x) || any(!nzchar(trimws(x))) || (unique && anyDuplicated(x)))
      stop(label," must contain nonempty literal text",if(unique) " without duplicates" else "",".")
  }
  exact_table(samples,c("Sample","Condition","Use","Experimental unit"),"Statistics samples")
  exact_table(comparisons,c("Name","Treatment","Control"),"Comparisons")
  use <- tolower(trimws(as.character(samples$Use)))
  if(anyNA(use) || any(!use %in% c("yes","no"))) stop("Use must be Yes or No for every sample.")
  selected <- samples[use=="yes",,drop=FALSE]
  if(!nrow(selected)) stop("Select independent experimental units with Use = Yes.")
  text_ids(selected$Sample,"Selected samples",TRUE)
  text_ids(selected$Condition,"Selected conditions")
  text_ids(selected[["Experimental unit"]],"Selected experimental units",TRUE)
  if(!nrow(comparisons)) stop("Declare at least one treatment/control comparison.")
  for(column in names(comparisons)) text_ids(comparisons[[column]],paste("Comparison",column),column=="Name")
  if(any(comparisons$Treatment==comparisons$Control)) stop("Treatment and Control must be different conditions.")
  if(anyDuplicated(comparisons[,c("Treatment","Control"),drop=FALSE])) stop("Repeated ordered condition comparisons are not permitted.")
  labels <- unique(c(comparisons$Treatment,comparisons$Control))
  if(any(!labels %in% selected$Condition)) stop("Comparisons refer to conditions without selected samples.")
  if(any(vapply(labels,function(label) sum(selected$Condition==label)<2L,logical(1))))
    stop("Each compared condition requires at least two declared independent experimental units.")

  static_observation <- c("observation_id","source_row","sequence","modifications","modifications_available",
    "candidate_entity_count","observation_excluded","observation_exclusion_reason")
  static_entity <- c("entity_id","entity_type","accessions","entity_source_row","entity_status","entity_excluded","entity_exclusion_reason")
  needed <- c(static_observation,static_entity,"sample_id","condition","log2_ratio","normalized_ratio","usable_pair","ratio_status","ratio_representation")
  if(!is.data.frame(ratios) || anyDuplicated(names(ratios)) || !all(needed %in% names(ratios)))
    stop("Statistics requires canonical per-sample ratio fields.")
  text_ids(ratios$observation_id,"Observation IDs");text_ids(ratios$sample_id,"Ratio sample IDs")
  if(!is.character(ratios$entity_id) || any(!is.na(ratios$entity_id) & !nzchar(trimws(ratios$entity_id))))
    stop("Entity IDs must be literal nonempty text or missing for absent links.")
  if(anyDuplicated(ratios[,c("observation_id","entity_id","sample_id"),drop=FALSE])) stop("Duplicate observation/entity/sample ratio rows.")
  if(!is.numeric(ratios$candidate_entity_count) || anyNA(ratios$candidate_entity_count) ||
      any(!is.finite(ratios$candidate_entity_count) | ratios$candidate_entity_count<0 | ratios$candidate_entity_count!=floor(ratios$candidate_entity_count)))
    stop("Candidate entity counts must be nonnegative integers.")
  if(!is.logical(ratios$usable_pair) || anyNA(ratios$usable_pair) || !is.numeric(ratios$log2_ratio) || !is.numeric(ratios$normalized_ratio))
    stop("Usable-pair flags and ratio quantities have invalid types.")
  if(!is.logical(ratios$observation_excluded) || anyNA(ratios$observation_excluded) || !is.logical(ratios$entity_excluded))
    stop("Canonical exclusion flags must be logical.")
  if(anyNA(ratios$entity_status) || any(!ratios$entity_status %in% c("matched","missing_protein_link","protein_entity_not_found")) ||
      any(is.na(ratios$entity_id)!=(ratios$entity_status=="missing_protein_link")) ||
      any((ratios$entity_status=="matched") & is.na(ratios$entity_excluded)))
    stop("Entity identities, states and exclusion flags are inconsistent.")
  text_ids(ratios$ratio_status,"Ratio statuses")
  valid_representation <- c("unavailable","linear_and_log2","log2_only_overflow","log2_only_underflow")
  if(anyNA(ratios$ratio_representation) || any(!ratios$ratio_representation %in% valid_representation) ||
      any(ratios$usable_pair != (ratios$ratio_status=="normalized")) ||
      any(ratios$usable_pair & (!is.finite(ratios$log2_ratio) | ratios$observation_excluded | ratios$entity_excluded %in% TRUE)) ||
      any(!ratios$usable_pair & !is.na(ratios$log2_ratio)) ||
      any(ratios$usable_pair == (ratios$ratio_representation=="unavailable")))
    stop("Ratio statuses, exclusions, representations and finite log2 values are inconsistent.")
  linear <- ratios$ratio_representation=="linear_and_log2"
  if(any(linear & (!is.finite(ratios$normalized_ratio) | ratios$normalized_ratio<=0)) || any(!linear & !is.na(ratios$normalized_ratio)))
    stop("Linear ratio values disagree with their declared representation.")
  if(any(ratios$usable_pair & (is.na(ratios$entity_id) | is.na(ratios$entity_excluded) | ratios$entity_status!="matched")))
    stop("Usable ratios require a matched, nonexcluded protein entity.")
  if(!is.character(ratios$condition) || anyNA(ratios$condition)) stop("Canonical ratio conditions must be text, using blank for unknown.")
  if(any(!selected$Sample %in% ratios$sample_id)) stop("Selected samples are absent from the ratio table.")
  source_samples <- unique(ratios$sample_id)
  for(sample in source_samples) {
    conditions <- unique(ratios$condition[ratios$sample_id==sample & nzchar(ratios$condition)])
    if(length(conditions)>1L) stop("Conflicting source condition labels for sample ",sample,".")
    position <- match(sample,selected$Sample)
    if(!is.na(position) && length(conditions) && conditions!=selected$Condition[position])
      stop("Statistics condition disagrees with the nonblank source condition for sample ",sample,".")
  }
  observation_index <- match(ratios$observation_id,unique(ratios$observation_id))
  entity_index <- match(ratios$entity_id,unique(ratios$entity_id))
  candidate_key <- paste(observation_index,entity_index,sep=":")
  keys <- unique(candidate_key)
  candidate_rows <- unname(split(seq_len(nrow(ratios)),factor(candidate_key,levels=keys)))
  same_fields <- function(rows,columns) all(vapply(columns,function(column) length(unique(ratios[[column]][rows]))==1L,logical(1)))
  for(rows in unname(split(seq_len(nrow(ratios)),observation_index))) {
    if(!same_fields(rows,static_observation)) stop("Observation metadata or declared candidate counts differ across ratio rows.")
    entities <- unique(ratios$entity_id[rows]);actual <- sum(!is.na(entities))
    if(actual!=ratios$candidate_entity_count[rows[1]] || (actual>0L && anyNA(entities)))
      stop("Declared candidate count disagrees with the retained entity links.")
  }
  for(rows in candidate_rows) {
    if(!same_fields(rows,static_entity)) stop("Protein-entity metadata differs across samples.")
    if(!setequal(ratios$sample_id[rows],source_samples)) stop("Canonical candidate rows must retain every source sample, including unavailable quantities.")
  }

  settings <- list(method="limma moderated two-group test of per-sample log2 peptide/protein ratios",design=design,
    response="log2_ratio",effect="mean treatment log2_ratio minus mean control log2_ratio",trend=FALSE,robust=FALSE,
    minimum_units_per_condition=2L,minimum_observations_for_moderation=2L,
    ambiguity="Only observations with exactly one declared protein-entity candidate enter the model and empirical Bayes pool.",
    variance="Exact zero and nonfinite within-condition residual variances are descriptive diagnostics; no variance floor is applied.",
    missingness="No imputation; at least two usable values in each condition are required for each tested observation.",
    multiple_testing="Benjamini-Hochberg correction within each contrast, over observations with valid tested p-values.",
    interpretation="Independent experimental units are user-declared; no paired, repeated, batch, site-level or protease-activity inference.",
    moderation="The minimum of two observations prevents a one-row unmoderated fallback; it does not establish prior stability or power.",
    infinite_prior="Infinite prior degrees of freedom are encoded as prior_df = NA and prior_df_infinite = TRUE; this is not a missing prior.",
    warnings=list(),contrast_families=list())
  if(!requireNamespace("limma",quietly=TRUE)) stop("Install limma for moderated ratio statistics.")
  settings$limma_version <- as.character(utils::packageVersion("limma"))
  result_list <- list();audit_list <- list()
  for(comparison_number in seq_len(nrow(comparisons))) {
    comparison <- comparisons[comparison_number,,drop=FALSE]
    metadata <- selected[selected$Condition %in% c(comparison$Control,comparison$Treatment),,drop=FALSE]
    is_control <- metadata$Condition==comparison$Control;is_treatment <- !is_control
    output <- ratios[vapply(candidate_rows,function(rows) rows[1],integer(1)),c(static_observation,static_entity),drop=FALSE]
    output$comparison <- comparison$Name;output$treatment <- comparison$Treatment;output$control <- comparison$Control
    for(column in c("n_control","n_treatment","eligible_observation_count","eb_pool_size","n_family_tests")) output[[column]] <- rep(0L,nrow(output))
    numeric_columns <- c("mean_control_log2_ratio","mean_treatment_log2_ratio","log2_effect","residual_df","residual_variance",
      "standard_error_unscaled","moderated_t","p_value","adjusted_p_value","posterior_variance","prior_variance","prior_df","total_df")
    for(column in numeric_columns) output[[column]] <- rep(NA_real_,nrow(output))
    output$prior_df_infinite <- rep(NA,nrow(output))
    output$eligible_for_fit <- rep(FALSE,nrow(output));output$tested <- rep(FALSE,nrow(output));output$status <- rep("insufficient_replicates",nrow(output))
    response <- matrix(NA_real_,nrow=nrow(output),ncol=nrow(metadata),dimnames=list(NULL,metadata$Sample))
    audits <- vector("list",nrow(output))
    for(i in seq_along(candidate_rows)) {
      rows <- candidate_rows[[i]]
      rows <- rows[match(metadata$Sample,ratios$sample_id[rows])]
      audit <- ratios[rows,,drop=FALSE]
      audit$comparison <- comparison$Name;audit$analysis_condition <- metadata$Condition
      audit$experimental_unit <- metadata[["Experimental unit"]]
      audit$used_for_fit <- FALSE;audit$analysis_status <- "insufficient_replicates"
      values <- audit$log2_ratio;values[!audit$usable_pair] <- NA_real_
      response[i,] <- values
      control_values <- values[is_control & is.finite(values)];treatment_values <- values[is_treatment & is.finite(values)]
      nc <- length(control_values);nt <- length(treatment_values)
      output$n_control[i] <- nc;output$n_treatment[i] <- nt
      if(nc) output$mean_control_log2_ratio[i] <- mean(control_values)
      if(nt) output$mean_treatment_log2_ratio[i] <- mean(treatment_values)
      if(nc && nt) output$log2_effect[i] <- mean(treatment_values)-mean(control_values)
      if(nc && nt && nc+nt>2L) {
        output$residual_df[i] <- nc+nt-2L
        output$residual_variance[i] <- (sum((control_values-mean(control_values))^2)+sum((treatment_values-mean(treatment_values))^2))/(nc+nt-2L)
      }
      status <- if(output$observation_excluded[i]) "excluded_observation" else if(output$candidate_entity_count[i]>1L) "ambiguous_baseline" else
        if(output$entity_status[i]!="matched") output$entity_status[i] else if(isTRUE(output$entity_excluded[i])) "excluded_protein_entity" else
        if(nc<2L || nt<2L) "insufficient_replicates" else
        if(!is.finite(output$log2_effect[i]) || !is.finite(output$residual_variance[i])) "numeric_range_exceeded" else
        if(output$residual_variance[i]==0) "zero_residual_variance" else "eligible"
      output$status[i] <- status;output$eligible_for_fit[i] <- status=="eligible"
      audits[[i]] <- audit
    }
    eligible <- which(output$eligible_for_fit)
    output$eligible_observation_count <- rep(length(eligible),nrow(output))
    if(length(eligible)<2L) {
      output$status[eligible] <- "insufficient_features_for_moderation"
    } else {
      design_matrix <- cbind(Control=as.numeric(is_control),Treatment=as.numeric(is_treatment))
      if(qr(design_matrix)$rank!=2L) stop("The comparison design is not full rank.")
      contrast <- matrix(c(-1,1),ncol=1L,dimnames=list(colnames(design_matrix),"Treatment-Control"))
      fit_warnings <- character()
      fit <- withCallingHandlers({
        fitted <- limma::lmFit(response[eligible,,drop=FALSE],design=design_matrix)
        fitted <- limma::contrasts.fit(fitted,contrasts=contrast)
        limma::eBayes(fitted,trend=FALSE,robust=FALSE)
      },warning=function(w) {fit_warnings <<- unique(c(fit_warnings,conditionMessage(w)));invokeRestart("muffleWarning")})
      settings$warnings[[comparison$Name]] <- fit_warnings
      output$eb_pool_size <- rep(length(eligible),nrow(output))
      output$residual_df[eligible] <- fit$df.residual
      output$residual_variance[eligible] <- fit$sigma^2
      output$standard_error_unscaled[eligible] <- fit$stdev.unscaled[,1]
      output$log2_effect[eligible] <- fit$coefficients[,1]
      output$moderated_t[eligible] <- fit$t[,1]
      output$p_value[eligible] <- fit$p.value[,1]
      output$posterior_variance[eligible] <- fit$s2.post
      output$prior_variance[eligible] <- rep(fit$s2.prior,length.out=length(eligible))
      output$prior_df[eligible] <- rep(fit$df.prior,length.out=length(eligible))
      output$prior_df_infinite[eligible] <- is.infinite(output$prior_df[eligible]) & output$prior_df[eligible]>0
      output$prior_df[eligible[output$prior_df_infinite[eligible]]] <- NA_real_
      output$total_df[eligible] <- fit$df.total
      valid <- is.finite(output$p_value[eligible]) & output$p_value[eligible]>=0 & output$p_value[eligible]<=1 &
        is.finite(output$moderated_t[eligible]) & is.finite(output$posterior_variance[eligible]) & output$posterior_variance[eligible]>0 &
        is.finite(output$total_df[eligible]) & output$total_df[eligible]>0
      output$tested[eligible] <- valid
      output$status[eligible] <- ifelse(valid,"tested","model_numeric_failure")
      output$p_value[eligible[!valid]] <- NA_real_
      tested <- eligible[valid]
      output$adjusted_p_value[tested] <- stats::p.adjust(output$p_value[tested],method="BH")
      output$n_family_tests <- rep(length(tested),nrow(output))
      for(i in eligible) audits[[i]]$used_for_fit <- audits[[i]]$usable_pair
    }
    for(i in seq_along(audits)) audits[[i]]$analysis_status <- output$status[i]
    for(column in numeric_columns) output[[column]][!is.finite(output[[column]])] <- NA_real_
    settings$contrast_families[[comparison$Name]] <- list(treatment=comparison$Treatment,control=comparison$Control,
      eligible_observation_count=length(eligible),eb_pool_size=if(length(eligible)>=2L) length(eligible) else 0L,
      n_family_tests=sum(output$tested))
    rownames(output) <- NULL
    result_list[[comparison_number]] <- output
    audit_list[[comparison_number]] <- do.call(rbind,audits)
  }
  results <- do.call(rbind,result_list);measurements <- do.call(rbind,audit_list)
  rownames(results) <- NULL;rownames(measurements) <- NULL
  list(results=results,measurements=measurements,samples=samples,comparisons=comparisons,settings=settings)
}
