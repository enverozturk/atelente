args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1]), winslash="/")))
source(file.path(root,"R","sample_metadata.R"),encoding="UTF-8")
source(file.path(root,"R","quantification.R"),encoding="UTF-8")
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))

fixture <- function(peptide = c(10,40), protein = c(20,40)) {
  columns <- paste0("Abundance: F",seq_along(peptide),": Sample")
  p <- data.frame(peptide_id="pep_1",source_row=2L,sequence="PEPTIDE",modifications="",
    mapping_status="one_match_in_scope",stringsAsFactors=FALSE)
  m <- data.frame(peptide_id="pep_1",accession="P1",is_master=TRUE,record_id=1L,start=2L,end=8L)
  measurements <- data.frame(peptide_id="pep_1",source_column=columns,abundance_type="abundance",
    measurement_label=paste0("F",seq_along(peptide)),value=as.numeric(peptide),original_value=as.character(peptide),
    value_status=ifelse(is.na(peptide),"missing","observed"),included_in_analysis=TRUE,stringsAsFactors=FALSE)
  proteins <- data.frame(Accession="P1","# Peptides"="7","# PSMs"="19",Contaminant=FALSE,check.names=FALSE)
  for(i in seq_along(columns)) proteins[[columns[i]]] <- as.character(protein[i])
  result <- list(peptides=p,matches=m,measurements=measurements,source_proteins=proteins)
  metadata <- pm_sample_template(result)
  metadata$samples$Condition <- rep(c("Control","Treatment"),each=length(peptide)/2)
  metadata$comparisons <- data.frame(Name="Treatment vs Control",Treatment="Treatment",Control="Control")
  list(result=result,metadata=metadata)
}

# One observation per condition is descriptive, with the expected multiplicative ratio.
x <- fixture(); before <- serialize(x,NULL)
out <- pm_adjust_peptide_effects(x$result,x$metadata)
stopifnot(identical(before,serialize(x,NULL)),nrow(out$effects)==1L,nrow(out$measurement_pairs)==2L,
  out$effects$status=="adjusted",out$effects$n_paired_treatment==1L,out$effects$n_paired_control==1L,
  out$effects$protein_source_row==2L,out$effects$protein_identified_peptide_count==7,
  out$effects$protein_psm_count==19,is.na(out$effects$protein_quantifying_peptide_count),
  !any(grepl("p.?value|confidence|significance",names(out$effects),ignore.case=TRUE)))
near(out$effects$raw_peptide_log2fc,2);near(out$effects$matched_peptide_log2fc,2)
near(out$effects$protein_log2fc,1);near(out$effects$adjusted_log2fc,log2((40/40)/(10/20)))
reverse <- x$metadata;reverse$comparisons$Treatment <- "Control";reverse$comparisons$Control <- "Treatment"
near(pm_adjust_peptide_effects(x$result,reverse)$effects$adjusted_log2fc,-1)
scaled <- fixture(c(100,400),c(2000,4000))
near(pm_adjust_peptide_effects(scaled$result,scaled$metadata)$effects$adjusted_log2fc,1)
fails(pm_adjust_peptide_effects(x$result,x$metadata,method="ratio_of_means"))

# A paired subset is shared by both terms, while the all-observation raw effect stays separate.
x4 <- fixture(c(4,64,16,1024),c(2,NA,8,NA));o4 <- pm_adjust_peptide_effects(x4$result,x4$metadata)
near(o4$effects$raw_peptide_log2fc,3);near(o4$effects$matched_peptide_log2fc,2)
near(o4$effects$protein_log2fc,2);near(o4$effects$adjusted_log2fc,0)
stopifnot(o4$effects$n_peptide_control==2,o4$effects$n_paired_control==1,
  o4$effects$n_selected_treatment==2,identical(o4$measurement_pairs$usable_pair,c(TRUE,FALSE,TRUE,FALSE)))

# Missing, zero, negative and invalid values remain distinct and never become pseudocounts.
for (case in list(list(value=NA_real_,state="missing"),list(value=0,state="zero"),
                 list(value=-1,state="negative"),list(value=Inf,state="invalid_numeric_value"))) {
  z <- fixture(c(10,case$value),c(20,40));oz <- pm_adjust_peptide_effects(z$result,z$metadata)
  stopifnot(is.na(oz$effects$adjusted_log2fc),oz$effects$status=="no_usable_treatment_pairs",
    oz$measurement_pairs$peptide_value_status[2]==case$state,!oz$measurement_pairs$usable_pair[2])
  if(identical(case$value,0)) stopifnot(oz$measurement_pairs$peptide_value[2]==0)
}
invalid <- fixture();invalid$result$source_proteins[["Abundance: F2: Sample"]] <- "not a number"
oi <- pm_adjust_peptide_effects(invalid$result,invalid$metadata)
stopifnot(oi$measurement_pairs$protein_value_status[2]=="invalid_numeric_value",is.na(oi$effects$adjusted_log2fc))
nonpositive <- fixture(c(10,40),c(0,-1));on <- pm_adjust_peptide_effects(nonpositive$result,nonpositive$metadata)
stopifnot(identical(on$measurement_pairs$protein_value_status,c("zero","negative")),
  on$effects$status=="no_usable_pairs_in_either_condition")
omitted <- fixture();omitted$result$measurements <- omitted$result$measurements[1,,drop=FALSE]
# Preserve column metadata via another peptide's retained measurement, as a missing row is not a missing cell.
omitted$result$peptides <- rbind(omitted$result$peptides,transform(omitted$result$peptides,peptide_id="pep_extra",mapping_status="no_exact_match"))
omitted$result$measurements <- rbind(omitted$result$measurements,transform(x$result$measurements[2,,drop=FALSE],peptide_id="pep_extra"))
oo <- pm_adjust_peptide_effects(omitted$result,omitted$metadata)
stopifnot(oo$measurement_pairs$peptide_value_status[2]=="measurement_not_exported",is.na(oo$effects$adjusted_log2fc))
excluded_measure <- fixture();excluded_measure$result$measurements$included_in_analysis[2] <- FALSE
oe <- pm_adjust_peptide_effects(excluded_measure$result,excluded_measure$metadata)
stopifnot(oe$measurement_pairs$peptide_value_status[2]=="excluded_measurement",is.na(oe$effects$raw_peptide_log2fc))

# Repeated sequence placements/FASTA records do not multiply the candidate baseline.
repeat_place <- fixture();repeat_place$result$matches <- rbind(repeat_place$result$matches,
  transform(repeat_place$result$matches,record_id=2L,start=12L,end=18L))
or <- pm_adjust_peptide_effects(repeat_place$result,repeat_place$metadata)
stopifnot(nrow(or$effects)==1L,nrow(or$measurement_pairs)==2L,or$effects$candidate_accession_count==1L)

# Multiple mapped masters are separate candidate cases. An exact isoform has its own baseline.
alternative <- repeat_place
alternative$result$matches <- rbind(alternative$result$matches,transform(alternative$result$matches[1,,drop=FALSE],accession="P1-2",record_id=3L))
isoform_protein <- alternative$result$source_proteins
isoform_protein$Accession <- "P1-2"
alternative$result$source_proteins <- rbind(alternative$result$source_proteins,isoform_protein)
alternative$result$source_proteins[["Abundance: F2: Sample"]][2] <- "80"
oa <- pm_adjust_peptide_effects(alternative$result,alternative$metadata)
stopifnot(nrow(oa$effects)==2L,nrow(oa$measurement_pairs)==4L,all(oa$effects$candidate_accession_count==2L))
near(oa$effects$adjusted_log2fc,c(1,0))
alternative$result$source_proteins <- alternative$result$source_proteins[1,,drop=FALSE]
missing_isoform <- pm_adjust_peptide_effects(alternative$result,alternative$metadata)
stopifnot(missing_isoform$effects$status[2]=="protein_accession_not_found",is.na(missing_isoform$effects$adjusted_log2fc[2]))

# Neither ambiguous duplicate protein rows nor non-master placements supply a guessed baseline.
duplicate <- fixture();duplicate$result$source_proteins <- rbind(duplicate$result$source_proteins,duplicate$result$source_proteins)
od <- pm_adjust_peptide_effects(duplicate$result,duplicate$metadata)
stopifnot(od$effects$status=="ambiguous_protein_rows",is.na(od$effects$adjusted_log2fc),is.na(od$effects$protein_source_row))
nonmaster <- fixture();nonmaster$result$matches$is_master <- FALSE
onm <- pm_adjust_peptide_effects(nonmaster$result,nonmaster$metadata)
stopifnot(onm$effects$status=="no_mapped_pd_master",onm$effects$candidate_accession_count==0L,is.na(onm$effects$accession))
flagged <- fixture();flagged$result$source_proteins$Contaminant <- "TRUE"
of <- pm_adjust_peptide_effects(flagged$result,flagged$metadata)
stopifnot(of$effects$status=="protein_flagged_contaminant",!any(of$measurement_pairs$usable_pair))

# Source observations, modifications and multiple comparisons remain independent keys.
multi <- fixture();multi$result$peptides <- rbind(multi$result$peptides,
  transform(multi$result$peptides,peptide_id="pep_2",source_row=3L,modifications="1xTest [P1]"))
multi$result$matches <- rbind(multi$result$matches,transform(multi$result$matches,peptide_id="pep_2"))
multi$result$measurements <- rbind(multi$result$measurements,transform(multi$result$measurements,peptide_id="pep_2",value=value*2))
multi$metadata$comparisons <- rbind(multi$metadata$comparisons,data.frame(Name="reverse",Treatment="Control",Control="Treatment"))
om <- pm_adjust_peptide_effects(multi$result,multi$metadata)
stopifnot(nrow(om$effects)==4L,nrow(om$measurement_pairs)==4L,length(unique(om$effects$peptide_id))==2L,
  all(om$effects$modifications[om$effects$peptide_id=="pep_2"]=="1xTest [P1]"))
near(om$effects$adjusted_log2fc,c(1,1,-1,-1))
multi$result$peptides$mapping_status[2] <- "excluded_contaminant"
om2 <- pm_adjust_peptide_effects(multi$result,multi$metadata)
stopifnot(nrow(om2$effects)==2L,!"pep_2" %in% om2$measurement_pairs$peptide_id)

# No retained observations returns schema-correct empty tables, while duplicate measurements fail.
empty <- fixture();empty$result$peptides$mapping_status <- "excluded_contaminant"
empty_out <- pm_adjust_peptide_effects(empty$result,empty$metadata)
stopifnot(identical(empty_out$effects,pm_empty_quant_effects()),identical(empty_out$measurement_pairs,pm_empty_quant_pairs()))
bad <- fixture();bad$result$measurements <- rbind(bad$result$measurements,bad$result$measurements[1,,drop=FALSE])
fails(pm_adjust_peptide_effects(bad$result,bad$metadata))
bad <- fixture();bad$result$source_proteins <- NULL;fails(pm_adjust_peptide_effects(bad$result,bad$metadata))
cat("Quantification tests passed.\n")
