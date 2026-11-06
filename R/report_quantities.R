# Use unique observation IDs so repeated placements and roles never multiply tests.
pm_paired_bh <- function(results, cleavage) {
  sites <- cleavage$sites
  paired <- sites$site_id[sites$neo_n_count > 0 & sites$neo_c_count > 0]
  evidence <- cleavage$evidence
  keep <- evidence$site_id %in% paired & (evidence$role %in% "spanning" |
    (evidence$role %in% c("N","C") & evidence$is_candidate %in% TRUE))
  results$paired_evidence <- results$peptide_id %in% evidence$peptide_id[keep]
  results$paired_adjusted_p_value <- rep(NA_real_, nrow(results))
  counts <- integer()
  for (contrast in unique(results$comparison)) {
    ix <- which(results$comparison == contrast & results$paired_evidence & results$tested)
    if (anyDuplicated(results$peptide_id[ix])) stop("Multiple tested baselines for one paired peptide in ", contrast)
    results$paired_adjusted_p_value[ix] <- stats::p.adjust(results$p_value[ix], method="BH")
    counts[contrast] <- length(ix)
  }
  list(results=results, counts=as.list(counts))
}

# Identity is established by the complete ordered source table, then source row,
# sequence and modification metadata. Sequence alone is never a join key.
pm_report_quantities <- function(result, ratios_dir, statistics_dir=NULL) {
  for(package in c("jsonlite","digest")) if(!requireNamespace(package,quietly=TRUE)) stop("Install ",package," first.")
  if(!exists("pm_read_table",mode="function")) stop("Source R/mapping.R before report attachment.")
  required <- function(x,fields,label) {
    if(!is.data.frame(x)||!all(fields %in% names(x))) stop("Missing ",label," fields: ",paste(setdiff(fields,names(x)),collapse=", "))
  }
  text <- function(x) { x<-as.character(x);x[!is.na(x)&x==""]<-NA_character_;x }
  same <- function(x,y) identical(text(x),text(y))
  same_table <- function(x,y) is.data.frame(x)&&is.data.frame(y)&&identical(names(x),names(y))&&nrow(x)==nrow(y)&&all(vapply(names(x),function(n)same(x[[n]],y[[n]]),logical(1)))
  key <- function(x,fields) {
    required(x,fields,"key")
    do.call(paste0,lapply(x[fields],function(v) {v<-text(v);ifelse(is.na(v),"N;",paste0("S",nchar(v,type="bytes"),":",v,";"))}))
  }
  unique_key <- function(x,fields,label) {if(anyDuplicated(key(x,fields))) stop("Duplicate ",label," keys.")}
  nonempty <- function(x,label) {if(anyNA(text(x))) stop(label," must be nonempty.")}
  numeric_fields <- c("source_row","entity_source_row","candidate_entity_count","peptide_input_value","protein_input_value",
    "peptide_log2_value","protein_log2_value","normalized_ratio","log2_ratio","n_control","n_treatment","eligible_observation_count",
    "eb_pool_size","n_family_tests","mean_control_log2_ratio","mean_treatment_log2_ratio","log2_effect","residual_df","residual_variance",
    "standard_error_unscaled","moderated_t","p_value","adjusted_p_value","posterior_variance","prior_variance","prior_df","total_df")
  logical_fields <- c("modifications_available","contaminant_flag","excluded","observation_excluded","entity_excluded","usable_pair",
    "prior_df_infinite","eligible_for_fit","tested","used_for_fit")
  typed <- function(x) {
    for(n in intersect(names(x),numeric_fields)) {
      raw<-x[[n]];v<-suppressWarnings(as.numeric(raw))
      if(any(!is.na(raw)&!is.finite(v))) stop("Invalid numeric report field: ",n)
      x[[n]]<-v
    }
    for(n in intersect(names(x),logical_fields)) {
      raw<-x[[n]]
      if(any(!is.na(raw)&!raw %in% c("TRUE","FALSE"))) stop("Invalid logical report field: ",n)
      x[[n]]<-ifelse(is.na(raw),NA,raw=="TRUE")
    }
    x
  }
  hash <- function(path) digest::digest(file=path,algo="sha256")
  read_manifest <- function(directory,name,stage) {
    path<-file.path(directory,name)
    if(!file.exists(path)) stop("Missing manifest: ",name)
    m<-jsonlite::read_json(path,simplifyVector=TRUE)
    if(!identical(m$stage,stage)) stop("Unexpected report input stage.")
    m
  }
  read_verified <- function(directory,manifest,name,convert=TRUE) {
    path<-file.path(directory,paste0(name,".csv"));expected<-manifest$csv_sha256[[name]]
    if(!file.exists(path)||is.null(expected)||!identical(expected,hash(path))) stop("Report input CSV hash mismatch: ",name)
    x<-pm_read_table(path)
    if(convert) typed(x) else x
  }
  required(result$peptides,c("peptide_id","source_row","sequence","modifications","modifications_available","mapping_status"),"mapping peptide")
  if(!is.data.frame(result$source_peptides)||nrow(result$source_peptides)!=nrow(result$peptides)) stop("Require original mapping source rows.")
  nonempty(result$peptides$peptide_id,"Mapping peptide IDs")
  unique_key(result$peptides,"peptide_id","mapping peptide");unique_key(result$peptides,"source_row","mapping source row")
  norm<-read_manifest(ratios_dir,"normalization_manifest.json","per_sample_ratios")
  source<-read_verified(ratios_dir,norm,"source_peptides",FALSE)
  # Normalization may consume the mapping source export with these two added
  # columns. They must be the actual mapping IDs and row positions, not labels
  # copied onto reordered data. No other columns are ignored.
  if(!same_table(source,result$source_peptides)) {
    prefix<-c("peptide_id","source_row")
    if(length(names(source))<2L||!identical(names(source)[1:2],prefix)||
       !same(source$peptide_id,result$peptides$peptide_id)||!same(source$source_row,result$peptides$source_row)||
       !same_table(source[,-(1:2),drop=FALSE],result$source_peptides))
      stop("Quantities and mapping must use the exact same ordered source peptide table.")
  }
  obs<-read_verified(ratios_dir,norm,"observations")
  required(obs,c("observation_id","source_row","sequence","modifications","modifications_available"),"quantity observation")
  nonempty(obs$observation_id,"Observation IDs")
  unique_key(obs,"observation_id","observation");unique_key(obs,"source_row","observation source row")
  if(nrow(obs)!=nrow(result$peptides)||anyNA(obs$source_row)||!setequal(obs$source_row,result$peptides$source_row)) stop("Observation source rows differ from mapping.")
  row_index<-match(obs$source_row,result$peptides$source_row)
  for(n in c("sequence","modifications","modifications_available"))
    if(!same(obs[[n]],result$peptides[[n]][row_index])) stop("Observation ",n," differs from mapping.")
  bridge<-data.frame(peptide_id=as.character(result$peptides$peptide_id[row_index]),obs[,c("observation_id","source_row","sequence","modifications","modifications_available"),drop=FALSE],
    mapping_status=as.character(result$peptides$mapping_status[row_index]),stringsAsFactors=FALSE)
  ratios<-read_verified(ratios_dir,norm,"sample_ratios")
  samples<-read_verified(ratios_dir,norm,"samples",FALSE)
  required(samples,c("sample_id","condition","peptide_column","protein_column"),"quantity sample")
  nonempty(samples$sample_id,"Sample IDs");unique_key(samples,"sample_id","sample")
  required(ratios,c("observation_id","entity_id","sample_id","source_row","sequence","modifications","modifications_available",
    "entity_type","accessions","candidate_entity_count","peptide_input_value","protein_input_value","peptide_log2_value","protein_log2_value",
    "normalized_ratio","log2_ratio","ratio_status","usable_pair","condition","peptide_column","protein_column",
    "peptide_scale","protein_scale","peptide_original_value","protein_original_value"),"ratio")
  nonempty(ratios$observation_id,"Ratio observation IDs");nonempty(ratios$sample_id,"Ratio sample IDs")
  unique_key(ratios,c("observation_id","entity_id","sample_id"),"ratio")
  sample_index<-match(ratios$sample_id,samples$sample_id)
  if(anyNA(sample_index)) stop("Ratio references an unknown sample.")
  for(n in c("condition","peptide_column","protein_column")) if(!same(ratios[[n]],samples[[n]][sample_index])) stop("Ratio sample metadata differs: ",n)
  attach <- function(table) {
    required(table,c("observation_id","source_row","sequence","modifications","modifications_available"),"attached observation")
    ix<-match(table$observation_id,bridge$observation_id)
    if(anyNA(ix)) stop("Unknown observation in attached results.")
    for(n in c("source_row","sequence","modifications","modifications_available")) if(!same(table[[n]],bridge[[n]][ix])) stop("Attached observation differs: ",n)
    table$peptide_id<-bridge$peptide_id[ix]
    table
  }
  ratios<-attach(ratios)
  if(!setequal(ratios$observation_id,bridge$observation_id)) stop("Missing observations in ratio table.")
  results<-measurements<-comparisons<-statistics_samples<-data.frame()
  statistics_settings<-NULL
  statistics_notes<-NULL
  provenance<-list(method="exact_ordered_source_table",normalization_manifest_sha256=hash(file.path(ratios_dir,"normalization_manifest.json")),
    ratio_sha256=norm$csv_sha256$sample_ratios,source_peptides_sha256=norm$csv_sha256$source_peptides)
  if(!is.null(statistics_dir)) {
    sm<-read_manifest(statistics_dir,"statistics_manifest.json","independent_unit_statistics")
    if(!identical(sm$inputs$ratios$sha256,provenance$ratio_sha256)||
       !identical(sm$inputs$ratio_manifest$sha256,provenance$normalization_manifest_sha256)) stop("Statistics were not generated from these normalized inputs.")
    results<-read_verified(statistics_dir,sm,"results")
    measurements<-read_verified(statistics_dir,sm,"measurements")
    comparisons<-read_verified(statistics_dir,sm,"comparisons",FALSE)
    statistics_samples<-read_verified(statistics_dir,sm,"samples",FALSE)
    required(comparisons,c("Name","Treatment","Control"),"comparison")
    nonempty(comparisons$Name,"Comparison names");unique_key(comparisons,"Name","comparison")
    required(statistics_samples,c("Sample","Condition","Use","Experimental unit"),"statistics sample")
    unique_key(statistics_samples,"Sample","statistics sample")
    required(results,c("comparison","entity_id","status","tested","treatment","control","log2_effect","p_value","adjusted_p_value","n_control","n_treatment"),"statistics result")
    required(measurements,c("comparison","entity_id","sample_id","used_for_fit","analysis_condition","experimental_unit","analysis_status"),"statistical measurement")
    unique_key(results,c("observation_id","entity_id","comparison"),"statistical result")
    unique_key(measurements,c("observation_id","entity_id","sample_id","comparison"),"statistical measurement")
    results<-attach(results);measurements<-attach(measurements)
    for(table in list(results,measurements)) if(any(!table$comparison %in% comparisons$Name)) stop("Unknown comparison in statistical table.")
    comparison_index<-match(results$comparison,comparisons$Name)
    if(!same(results$treatment,comparisons$Treatment[comparison_index])||!same(results$control,comparisons$Control[comparison_index])) stop("Statistical contrast labels disagree.")
    selected<-statistics_samples[tolower(trimws(statistics_samples$Use)) %in% "yes",,drop=FALSE]
    selected_index<-match(measurements$sample_id,selected$Sample)
    if(anyNA(selected_index)||!same(measurements$analysis_condition,selected$Condition[selected_index])||
       !same(measurements$experimental_unit,selected[["Experimental unit"]][selected_index])) stop("Statistical sample metadata disagrees with measurements.")
    comparison_index<-match(measurements$comparison,comparisons$Name)
    if(anyNA(measurements$analysis_condition)||any(!((measurements$analysis_condition==comparisons$Treatment[comparison_index])|
       (measurements$analysis_condition==comparisons$Control[comparison_index])))) stop("Statistical sample is outside its compared conditions.")
    # Statistics source values must remain the same exported ratios, including
    # NA/zero distinctions; only analysis fields are additional.
    rkey<-key(ratios,c("observation_id","entity_id","sample_id"))
    ix<-match(key(measurements,c("observation_id","entity_id","sample_id")),rkey)
    if(anyNA(ix)||!all(names(ratios) %in% names(measurements))) stop("Statistical measurements do not match ratio identities.")
    for(n in names(ratios)) if(!same(measurements[[n]],ratios[[n]][ix])) stop("Statistical source value differs: ",n)
    candidates<-unique(ratios[,c("observation_id","entity_id"),drop=FALSE])
    if(any(!key(results,c("observation_id","entity_id")) %in% key(candidates,c("observation_id","entity_id")))) stop("Statistics references an unknown quantitative candidate.")
    if(nrow(results)!=nrow(candidates)*nrow(comparisons)) stop("Statistical candidate results are incomplete.")
    expected_measurements<-nrow(candidates)*sum(vapply(seq_len(nrow(comparisons)),function(i)
      sum(selected$Condition %in% c(comparisons$Treatment[i],comparisons$Control[i])),integer(1)))
    if(nrow(measurements)!=expected_measurements) stop("Statistical sample measurements are incomplete.")
    for(n in c("p_value","adjusted_p_value")) if(any(!is.na(results[[n]])&(results[[n]]<0|results[[n]]>1))) stop("Invalid statistical probability: ",n)
    if(anyNA(results$tested)||anyNA(results$status)||any(results$tested != (results$status=="tested"))||
       any(results$tested & (is.na(results$p_value)|is.na(results$adjusted_p_value)))||
       any(!results$tested & (!is.na(results$p_value)|!is.na(results$adjusted_p_value)))) stop("Statistics test status and probabilities disagree.")
    ri<-match(key(measurements,c("observation_id","entity_id","comparison")),key(results,c("observation_id","entity_id","comparison")))
    if(anyNA(ri)||!same(measurements$analysis_status,results$status[ri])) stop("Statistical sample status disagrees with result.")
    statistics_settings<-sm$settings
    statistics_notes<-sm$notes
    provenance$statistics_manifest_sha256<-hash(file.path(statistics_dir,"statistics_manifest.json"))
  }
  provenance$row_counts<-list(observations=nrow(bridge),ratios=nrow(ratios),results=nrow(results),measurements=nrow(measurements))
  list(observations=bridge,ratios=ratios,samples=samples,results=results,measurements=measurements,
    comparisons=comparisons,statistics_samples=statistics_samples,
    settings=list(quantity=norm$settings,spec=norm$spec,statistics=statistics_settings,notes=list(quantity=norm$notes,statistics=statistics_notes)),provenance=provenance)
}
