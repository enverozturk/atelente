# Workflow adapter for the existing experimental measurement-pair model.
# Mapping links are annotations of measurements, never additional independent tests.
pm_pair_workflow <- function(result, design, representative_column=NULL,
                             assignment_mode=c("unique","representative"), minimum=3L,
                             test=TRUE, experimental_units=NULL, abundance_type=NULL,
                             n_terminal_states=NULL) {
  assignment_mode <- match.arg(assignment_mode)
  if(!is.logical(test)||length(test)!=1L||is.na(test)) stop("test must be TRUE or FALSE.")
  if(!is.numeric(minimum)||length(minimum)!=1L||is.na(minimum)||minimum<3L||minimum!=as.integer(minimum))
    stop("Pair inference minimum must be an integer >=3.")
  if(is.null(result$cleavage)) stop("Pair workflow requires cleavage evidence.")
  if(!is.list(design)||!is.data.frame(design$samples)||!is.data.frame(design$comparisons))
    stop("Explicit samples and comparisons are required.")
  fields <- c("Sample","Condition","Use","Peptide column")
  if(!all(fields %in% names(design$samples))||!all(c("Name","Treatment","Control") %in% names(design$comparisons)))
    stop("Pair workflow sample/comparison columns are missing.")
  if(anyNA(design$samples$Use)||any(!design$samples$Use %in% c("Yes","No")))
    stop("Sample Use must explicitly be Yes or No.")
  samples <- design$samples[design$samples$Use=="Yes",,drop=FALSE]
  for(field in fields) if(anyNA(samples[[field]])||any(!nzchar(as.character(samples[[field]]))))
    stop("Included sample fields must be nonempty.")
  if(!nrow(samples)||anyDuplicated(samples$Sample)||anyDuplicated(samples[["Peptide column"]]))
    stop("Include unique samples and unique quantity columns.")
  comparisons <- design$comparisons
  for(field in c("Name","Treatment","Control")) if(anyNA(comparisons[[field]])||any(!nzchar(as.character(comparisons[[field]]))))
    stop("Comparison fields must be nonempty.")
  if(!nrow(comparisons)||anyDuplicated(comparisons$Name)||
     anyDuplicated(comparisons[,c("Treatment","Control"),drop=FALSE])||any(comparisons$Treatment==comparisons$Control)||
     !all(c(comparisons$Treatment,comparisons$Control) %in% samples$Condition))
    stop("Declare unique comparisons between distinct included conditions.")
  p <- result$peptides;m <- result$matches;source <- result$source_peptides
  if(nrow(p)!=nrow(source)||anyDuplicated(p$peptide_id)) stop("Source observation identity is inconsistent.")
  ids <- p$peptide_id
  if(assignment_mode=="representative") {
    if(!is.character(representative_column)||length(representative_column)!=1L||is.na(representative_column)||
       !representative_column %in% names(source)) stop("Choose an existing representative-protein column.")
    representatives <- setNames(lapply(as.character(source[[representative_column]]),function(x) sort(unique(pm_split(x)))),ids)
  } else representatives <- NULL
  candidate_ids <- split(as.character(result$candidates$accession),result$candidates$peptide_id)
  mapped <- split(seq_len(nrow(m)),m$peptide_id)
  # Infer only a named, exact N-terminal modification. Unknown chemistry remains
  # its own exploratory class; the PKD adapter may supply its verified form audit.
  states <- setNames(rep("unknown",length(ids)),ids)
  label <- result$cleavage$profile$pre_digest_n_label
  for(i in seq_along(ids)) {
    if(!isTRUE(p$modifications_available[i])||is.na(p$modifications[i])) next
    mod <- p$modifications[i]
    if(!is.null(label)&&isTRUE(pm_n_label(mod,TRUE,label))) states[i] <- "labelled"
    else if(isTRUE(pm_n_label(mod,TRUE,"Acetyl"))) states[i] <- "acetyl"
    else if(!grepl("N-Term",mod,fixed=TRUE)) states[i] <- "unlabelled"
  }
  if(!is.null(n_terminal_states)) {
    if(!is.data.frame(n_terminal_states)||!all(c("peptide_id","state") %in% names(n_terminal_states))||
       anyDuplicated(n_terminal_states$peptide_id)||!all(n_terminal_states$peptide_id %in% ids))
      stop("N-terminal state audit requires unique known peptide_id and state columns.")
    state <- as.character(n_terminal_states$state)
    state[state=="TMT"] <- "labelled";state[state=="Acetyl"] <- "acetyl"
    if(anyNA(state)||any(!state %in% c("labelled","acetyl","unlabelled","unknown")))
      stop("Unknown N-terminal state in supplied audit.")
    # A verified audit may resolve unavailable or unnamed chemistry, but cannot
    # overwrite a recognized source state or invent a label on a known blank form.
    at <- match(n_terminal_states$peptide_id,ids)
    if(any(states[at]!="unknown" & states[at]!=state))
      stop("N-terminal state audit contradicts a known source state.")
    states[at] <- state
  }
  if(length(setdiff(samples[["Peptide column"]],unique(result$measurements$source_column))))
    stop("Unknown selected quantity column.")
  selected_type <- unique(result$measurements$abundance_type[result$measurements$source_column %in% samples[["Peptide column"]]])
  if(length(selected_type)!=1L||!selected_type %in% c("abundance","normalized"))
    stop("Pair workflow requires individual sample columns of one abundance type.")
  if(!is.null(abundance_type)&&!identical(abundance_type,selected_type))
    stop("abundance_type conflicts with selected individual sample columns.")
  abundance_type <- selected_type
  ratios <- pm_neo_spanning_ratios(result$cleavage$evidence,result$measurements,abundance_type)
  ratios <- ratios[ratios$source_column %in% samples[["Peptide column"]],,drop=FALSE]
  contrasts <- pm_neo_spanning_contrasts(ratios,design)
  keys <- c("site_id","record_id","accession","bond_left","bond_right","role","neo_peptide_id","spanning_peptide_id")
  links <- unique(ratios[,keys,drop=FALSE]);rownames(links) <- NULL
  pair_key <- function(neo,span) paste(neo,span,sep="/")
  directed <- pair_key(links$neo_peptide_id,links$spanning_peptide_id)
  distinct <- unique(directed)
  hash <- vapply(distinct,function(key) paste0("pair_",digest::digest(key,algo="sha256",serialize=FALSE)),character(1))
  if(anyDuplicated(hash)) stop("Measurement pair hash collision.")
  pair_ids <- setNames(hash,distinct)
  links$measurement_pair_id <- unname(pair_ids[directed])
  links$neo_nt_state <- unname(states[links$neo_peptide_id])
  n_family <- if(is.null(label)) "N_standard" else if(label %in% c("TMT","TMTpro")) "N_TMT_primary" else "N_labelled_primary"
  nf <- if(is.null(label)) rep("N_standard",nrow(links)) else
    ifelse(links$neo_nt_state=="labelled",n_family,
      ifelse(links$neo_nt_state=="acetyl","N_acetyl_exploratory",
        ifelse(links$neo_nt_state=="unlabelled","N_unlabelled_exploratory","N_unresolved_exploratory")))
  links$family <- ifelse(links$role=="C","C_secondary",nf)
  # The vector branch above deliberately distinguishes absent and unknown labels.
  if(!nrow(links)) links$family <- character()
  links$assignment_mode <- rep(assignment_mode,nrow(links))
  links$representative_accession <- rep("",nrow(links))
  links$common_master_accessions <- rep("",nrow(links))
  links$common_master_count <- integer(nrow(links))
  links$placement_policy <- rep("first_sequence_occurrence",nrow(links))
  links$neo_placement_count <- integer(nrow(links))
  links$spanning_placement_count <- integer(nrow(links))
  for(field in c("neo_selected_start","neo_selected_end","spanning_selected_start","spanning_selected_end"))
    links[[field]] <- rep(NA_integer_,nrow(links))
  links$assignment_eligible <- rep(FALSE,nrow(links))
  links$scope_eligible <- rep(TRUE,nrow(links));links$exclusion_reason <- rep("",nrow(links))
  evidence <- result$cleavage$evidence
  evidence <- evidence[evidence$role=="spanning" | evidence$is_candidate %in% TRUE,,drop=FALSE]
  evidence_key <- function(site,peptide,role,start,end) paste(site,peptide,role,start,end,sep="/")
  evidence_keys <- unique(evidence_key(evidence$site_id,evidence$peptide_id,evidence$role,evidence$start,evidence$end))
  evidence_lookup <- list2env(as.list(setNames(rep(TRUE,length(evidence_keys)),evidence_keys)),parent=emptyenv())
  # Same per-link rules as before, on plain vectors (no per-row data.frame
  # indexing); results are written back to `links` once after the loop.
  nl <- nrow(links)
  l_site <- links$site_id; l_neo <- links$neo_peptide_id; l_span <- links$spanning_peptide_id
  l_acc <- links$accession; l_rec <- links$record_id; l_role <- links$role
  l_left <- links$bond_left; l_right <- links$bond_right
  m_acc <- m$accession; m_rec <- m$record_id; m_start <- m$start; m_end <- m$end
  r_common <- links$common_master_accessions; r_ncommon <- links$common_master_count
  r_rep <- links$representative_accession
  r_nneo <- links$neo_placement_count; r_nspan <- links$spanning_placement_count
  r_ns <- links$neo_selected_start; r_ne <- links$neo_selected_end
  r_ss <- links$spanning_selected_start; r_se <- links$spanning_selected_end
  r_assign <- links$assignment_eligible; r_scope <- links$scope_eligible; r_reason <- links$exclusion_reason
  l_ref <- result$fasta$sequence[match(l_rec,result$fasta$record_id)]
  all_scope <- identical(result$settings$scope,"all")
  # Unique (record,start,end) placements of one form on one accession, first-seen order.
  placements_on <- function(id,acc) {
    ii <- mapped[[id]]; ii <- ii[m_acc[ii]==acc]
    ii[!duplicated(paste(m_rec[ii],m_start[ii],m_end[ii],sep="\r"))]
  }
  for(i in seq_len(nl)) {
    neo <- l_neo[i];span <- l_span[i]
    reason <- character();acc <- l_acc[i]
    if(assignment_mode=="representative") {
      a <- representatives[[neo]];b <- representatives[[span]]
      common <- sort(intersect(a,b))
      r_common[i] <- paste(common,collapse=";")
      r_ncommon[i] <- length(common)
      if(!length(a)||!length(b)) reason <- "Missing representative assignment"
      else if(!length(common)) reason <- "No shared master accession"
      else if(!acc %in% common) reason <- "Outside shared master assignments"
      else r_rep[i] <- acc
    } else {
      a <- unique(candidate_ids[[neo]]);b <- unique(candidate_ids[[span]])
      # All-proteome mapping has no declared candidate list: use all exact hits.
      if(all_scope) {
        a <- unique(m_acc[mapped[[neo]]]);b <- unique(m_acc[mapped[[span]]])
      }
      if(length(a)!=1L||length(b)!=1L||!identical(a,b)||!identical(a,acc))
        reason <- "Ambiguous protein assignment"
      else r_rep[i] <- acc
    }
    if(!length(reason)) {
      # Match TermineR's first sequence occurrence convention within one record.
      # Never choose the first accession or silently resolve duplicate records.
      ni <- placements_on(neo,acc); si <- placements_on(span,acc)
      r_nneo[i] <- length(ni); r_nspan[i] <- length(si)
      if(!length(ni)||!length(si)||length(unique(m_rec[ni]))!=1L||length(unique(m_rec[si]))!=1L||
         m_rec[ni[1]]!=l_rec[i]||m_rec[si[1]]!=l_rec[i])
        reason <- "Ambiguous placement on representative"
      else {
        n1 <- ni[order(m_start[ni],m_end[ni])][1]; s1 <- si[order(m_start[si],m_end[si])][1]
        r_ns[i] <- m_start[n1];r_ne[i] <- m_end[n1]
        r_ss[i] <- m_start[s1];r_se[i] <- m_end[s1]
        valid_neo <- if(l_role[i]=="N") m_start[n1]==l_right[i] else m_end[n1]==l_left[i]
        valid_span <- m_start[s1]<=l_left[i] & m_end[s1]>=l_right[i]
        if(!any(valid_neo)||!any(valid_span)) {
          reason <- "First sequence occurrences do not support this bond"
        } else if(!exists(evidence_key(l_site[i],neo,l_role[i],m_start[n1],m_end[n1]),envir=evidence_lookup,inherits=FALSE)||
                  !exists(evidence_key(l_site[i],span,"spanning",m_start[s1],m_end[s1]),envir=evidence_lookup,inherits=FALSE)) {
          reason <- "First sequence occurrence is outside candidate evidence scope"
        }
      }
    }
    r_assign[i] <- !length(reason)
    ref <- l_ref[i]
    if(l_left[i]==1L && !is.na(ref) && startsWith(ref,"M"))
      reason <- c(reason,"Possible initiator-Met processing; outside internal test scope")
    if(length(reason)) {r_scope[i] <- FALSE;r_reason[i] <- paste(reason,collapse="; ")}
  }
  links$common_master_accessions <- r_common; links$common_master_count <- r_ncommon
  links$representative_accession <- r_rep
  links$neo_placement_count <- r_nneo; links$spanning_placement_count <- r_nspan
  links$neo_selected_start <- r_ns; links$neo_selected_end <- r_ne
  links$spanning_selected_start <- r_ss; links$spanning_selected_end <- r_se
  links$assignment_eligible <- r_assign; links$scope_eligible <- r_scope; links$exclusion_reason <- r_reason
  # Multiple valid attributions remain links to one measured pair, never tests.
  group <- paste(links$measurement_pair_id,links$family,sep="/")
  eligible_counts <- table(group[links$assignment_eligible])
  links$eligible_mapping_count <- as.integer(eligible_counts[group])
  links$eligible_mapping_count[is.na(links$eligible_mapping_count)] <- 0L
  links$attribution_ambiguous <- links$common_master_count>1L | links$eligible_mapping_count>1L |
    links$neo_placement_count>1L | links$spanning_placement_count>1L |
    grepl("Ambiguous placement",links$exclusion_reason,fixed=TRUE)
  # Verify that re-used mappings carry exactly one value per pair and sample.
  rid <- unname(pair_ids[pair_key(ratios$neo_peptide_id,ratios$spanning_peptide_id)])
  sample_index <- match(ratios$source_column,samples[["Peptide column"]])
  measurement <- paste(rid,sample_index,sep="/")
  first <- !duplicated(measurement)
  first_index <- match(measurement,measurement[first])
  a <- ratios$log2_ratio;b <- a[first][first_index]
  if(any(xor(is.na(a),is.na(b)))||any(a[!is.na(a)]!=b[!is.na(b)]))
    stop("Alternative mappings contain inconsistent measurement-pair values.")
  response <- matrix(NA_real_,length(distinct),nrow(samples),dimnames=list(unname(pair_ids),samples$Sample))
  if(length(rid)) response[cbind(match(rid[first],rownames(response)),sample_index[first])] <- a[first]
  membership <- unique(links[links$scope_eligible,c("measurement_pair_id","family"),drop=FALSE])
  metadata <- data.frame(sample=as.character(samples$Sample),condition=as.character(samples$Condition),stringsAsFactors=FALSE)
  if(test) {
    units <- if(is.null(experimental_units)&&"Experimental unit" %in% names(samples)) samples[,c("Sample","Experimental unit"),drop=FALSE] else experimental_units
    if(!is.data.frame(units)||!all(c("Sample","Experimental unit") %in% names(units))||
       anyDuplicated(units$Sample)||!all(samples$Sample %in% units$Sample))
      stop("Pair tests require explicit independent experimental units for every included sample.")
    unit_rows <- match(samples$Sample,units$Sample)
    for(field in intersect(c("Condition","Use"),names(units))) {
      supplied <- as.character(units[[field]][unit_rows])
      if(anyNA(supplied)||!identical(supplied,as.character(samples[[field]])))
        stop("Experimental-unit ",field," conflicts with included sample design.")
    }
    metadata$experimental_unit <- as.character(units[["Experimental unit"]][match(samples$Sample,units$Sample)])
    if(anyNA(metadata$experimental_unit)||any(!nzchar(metadata$experimental_unit))||anyDuplicated(metadata$experimental_unit))
      stop("Pair tests require distinct nonempty independent experimental units.")
  }
  fits <- list();reports <- list();exclusions <- list();warnings <- list()
  for(i in seq_len(nrow(comparisons))) {
    cmp <- comparisons[i,,drop=FALSE]
    ex <- links[!links$scope_eligible,,drop=FALSE];ex$comparison <- rep(cmp$Name,nrow(ex));ex$reason <- ex$exclusion_reason
    exclusions[[i]] <- ex
    if(test&&nrow(membership)) {
      fit <- pm_test_pair_matrix(response,metadata,membership,cmp$Control,cmp$Treatment,minimum)
      pm_check_pair_arithmetic(fit,response,metadata,cmp$Control,cmp$Treatment)
      fit$results$comparison <- cmp$Name
      fits[[i]] <- fit$results;warnings[[cmp$Name]] <- fit$warnings
      selected <- links[links$scope_eligible,,drop=FALSE]
      selected$comparison <- cmp$Name
      reports[[i]] <- merge(selected,fit$results,by=c("measurement_pair_id","family","comparison"),sort=FALSE)
    }
  }
  bind <- function(x) if(length(x)) {y<-do.call(rbind,x);rownames(y)<-NULL;y} else data.frame()
  list(ratios=ratios,contrasts=contrasts,links=links,excluded=bind(exclusions),
    measurement_matrix=response,membership=membership,results=bind(reports),measurement_results=bind(fits),
    samples=metadata,warnings=warnings,
    settings=list(scope_version="first_occurrence_pair_scope_v5",analysis_scope=result$cleavage$profile$analysis_scope,
      assignment_mode=assignment_mode,representative_column=representative_column,minimum=minimum,test=test,
      placement_policy="First exact sequence occurrence per form within a single accession/FASTA record; repeated positions retained in audit. This is a convention, not resolved localization.",
      natural_start_policy="Initial-M bond 1|2 excluded from internal inference; descriptive results retained.",
      unit="Distinct directed observation pair once per family and contrast; protein links are annotations.",
      method=if(test) "OLS limma; BH within family and contrast" else "Descriptive only"))
}
