args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1]), winslash="/")))
source(file.path(root,"R","sample_ratios.R"),encoding="UTF-8")
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))

# Independent canonical fixtures: arithmetic tests do not reuse the importer.
fixture <- function(peptide=c(10,40),protein=c(20,40),peptide_scale="linear",protein_scale="linear") {
  stopifnot(length(peptide)==length(protein))
  samples <- data.frame(sample_id=paste0("S",seq_along(peptide)),condition="Control",
    peptide_column=paste0("Peptide ",seq_along(peptide)),protein_column=paste0("Protein ",seq_along(peptide)),stringsAsFactors=FALSE)
  measurements <- function(values,scale,layer) {
    state <- ifelse(is.na(values),"missing",ifelse(!is.finite(values),"invalid_numeric_value",
      ifelse(scale=="linear" & values==0,"zero",ifelse(scale=="linear" & values<0,"negative","observed"))))
    logvalue <- rep(NA_real_,length(values));valid <- state=="observed"
    logvalue[valid] <- if(scale=="linear") log2(values[valid]) else values[valid]
    out <- data.frame(row_id=if(layer=="peptide") "o1" else "e1",sample_id=samples$sample_id,
      source_column=samples[[paste0(layer,"_column")]],input_value=as.numeric(values),original_value=as.character(values),
      value_status=state,log2_value=logvalue,stringsAsFactors=FALSE)
    names(out)[1] <- if(layer=="peptide") "observation_id" else "entity_id"
    out
  }
  list(observations=data.frame(observation_id="o1",source_row=2L,sequence="PEPTIDE",modifications="",
      modifications_available=TRUE,excluded=FALSE,exclusion_reason="",stringsAsFactors=FALSE),
    entities=data.frame(entity_id="e1",entity_type="protein",accessions="P1",source_row=7L,
      excluded=FALSE,exclusion_reason="",stringsAsFactors=FALSE),
    links=data.frame(observation_id="o1",entity_id="e1",stringsAsFactors=FALSE),samples=samples,
    peptide_measurements=measurements(peptide,peptide_scale,"peptide"),protein_measurements=measurements(protein,protein_scale,"protein"),
    spec=list(peptide=list(scale=peptide_scale,quantity_type="intensity"),protein=list(scale=protein_scale,quantity_type="intensity")),
    source_peptides=data.frame(sequence="PEPTIDE",note="unchanged"),source_proteins=data.frame(entity="P1"))
}

# Same-condition replicas need no contrast and retain source quantities and rows.
x <- fixture();before <- serialize(x,NULL);out <- pm_sample_ratios(x);r <- out$ratios
stopifnot(identical(before,serialize(x,NULL)),nrow(r)==2L,all(r$condition=="Control"),all(r$usable_pair),
  all(r$ratio_status=="normalized"),all(r$ratio_representation=="linear_and_log2"),
  all(r$source_row==2L),all(r$entity_source_row==7L),all(r$candidate_entity_count==1L),
  all(r$modifications_available),out$summary$peptide_measurements==2L,
  !any(grepl("p.?value|fold.?change|significance",names(r),ignore.case=TRUE)))
near(r$normalized_ratio,c(.5,1));near(r$log2_ratio,c(-1,0))
stopifnot(identical(r$peptide_original_value,c("10","40")),identical(r$protein_original_value,c("20","40")))
x$samples$condition <- "";stopifnot(nrow(pm_sample_ratios(x)$ratios)==2L)

# Expected 1/2/4 quantities on linear, log2 and mixed declared scales.
linear <- pm_sample_ratios(fixture(c(10,40,120),c(10,20,30)))$ratios
logged <- pm_sample_ratios(fixture(log2(c(10,40,120)),log2(c(10,20,30)),"log2","log2"))$ratios
mixed <- pm_sample_ratios(fixture(log2(c(10,40,120)),c(10,20,30),"log2","linear"))$ratios
near(linear$normalized_ratio,c(1,2,4));near(linear$log2_ratio,c(0,1,2))
near(logged$normalized_ratio,linear$normalized_ratio);near(mixed$normalized_ratio,linear$normalized_ratio)
r <- pm_sample_ratios(fixture(c(0,-1,-2),c(-1,0,-3),"log2","log2"))$ratios
stopifnot(all(r$usable_pair),all(r$peptide_value_status=="observed"));near(r$normalized_ratio,c(2,.5,2))

# No zero substitution, pseudocount, or conversion of bad data to observations.
r <- pm_sample_ratios(fixture(c(NA,0,-2,Inf,1,1,1,1),c(1,1,1,1,NA,0,-2,Inf)))$ratios
stopifnot(identical(r$ratio_status,c("peptide_missing","peptide_zero","peptide_negative","peptide_invalid_numeric_value",
  "protein_missing","protein_zero","protein_negative","protein_invalid_numeric_value")),
  all(!r$usable_pair),all(is.na(r$normalized_ratio)),all(is.na(r$log2_ratio)))
x <- fixture();x$protein_measurements <- x$protein_measurements[-1,,drop=FALSE]
stopifnot(pm_sample_ratios(x)$ratios$ratio_status[1]=="protein_measurement_not_exported")
x <- fixture();x$peptide_measurements <- x$peptide_measurements[-1,,drop=FALSE]
stopifnot(pm_sample_ratios(x)$ratios$ratio_status[1]=="peptide_measurement_not_exported")

# Preserve finite log2 results even outside the linear representable range.
r <- pm_sample_ratios(fixture(c(2000,-2000,.Machine$double.xmax),c(0,0,-.Machine$double.xmax),"log2","log2"))$ratios
stopifnot(identical(r$usable_pair,c(TRUE,TRUE,FALSE)),all(is.na(r$normalized_ratio)),
  identical(r$ratio_representation,c("log2_only_overflow","log2_only_underflow","unavailable")),
  identical(r$ratio_status,c("normalized","normalized","numeric_range_exceeded")))
near(r$log2_ratio,c(2000,-2000,NA_real_))
r <- pm_sample_ratios(fixture(.Machine$double.xmax,.Machine$double.xmin))$ratios
stopifnot(r$usable_pair,r$ratio_representation=="log2_only_overflow",is.finite(r$log2_ratio))

# A group is one denominator, regardless of how many accessions it carries.
x <- fixture();x$entities$entity_type <- "protein_group";x$entities$accessions <- "P1;P1-2;P2"
r <- pm_sample_ratios(x)$ratios
stopifnot(nrow(r)==2L,all(r$entity_type=="protein_group"),all(r$accessions=="P1;P1-2;P2"))
# Alternative candidates reuse the same observation, not new quantitative evidence.
e2 <- x$entities;e2$entity_id <- "e2";e2$accessions <- "P3";e2$source_row <- 8L
x$entities <- rbind(x$entities,e2)
m2 <- x$protein_measurements;m2$entity_id <- "e2";m2$input_value <- m2$input_value*2;m2$log2_value <- m2$log2_value+1
x$protein_measurements <- rbind(x$protein_measurements,m2)
x$links <- rbind(x$links,data.frame(observation_id="o1",entity_id="e2"))
out <- pm_sample_ratios(x);r <- out$ratios
stopifnot(nrow(r)==4L,all(r$candidate_entity_count==2L),length(unique(r$observation_id))==1L,
  out$summary$observations==1L,out$summary$peptide_measurements==2L)
near(r$normalized_ratio,c(.5,1,.25,.5))

# Missing references are visible diagnostics; explicit exclusions never yield ratios.
x <- fixture();x$links <- x$links[FALSE,,drop=FALSE];r <- pm_sample_ratios(x)$ratios
stopifnot(nrow(r)==2L,all(r$ratio_status=="missing_protein_link"),all(r$candidate_entity_count==0L),all(is.na(r$entity_id)))
x <- fixture();x$links$entity_id <- "unknown";r <- pm_sample_ratios(x)$ratios
stopifnot(nrow(r)==2L,all(r$ratio_status=="protein_entity_not_found"),all(r$entity_id=="unknown"),all(!r$usable_pair))
x <- fixture();x$entities$excluded <- TRUE;x$entities$exclusion_reason <- "source_contaminant_flag"
r <- pm_sample_ratios(x)$ratios
stopifnot(all(r$ratio_status=="excluded_protein_entity"),all(!r$usable_pair),all(r$entity_exclusion_reason=="source_contaminant_flag"))
x$observations$excluded <- TRUE;x$observations$exclusion_reason <- "linked_contaminant_entity"
r <- pm_sample_ratios(x)$ratios
stopifnot(all(r$ratio_status=="excluded_observation"),all(!r$usable_pair),all(r$observation_exclusion_reason=="linked_contaminant_entity"))

# Distinct modification observations with the same sequence remain distinct.
x <- fixture();o2 <- x$observations;o2$observation_id <- "o2";o2$source_row <- 3L;o2$modifications <- "Oxidation [M1]"
x$observations <- rbind(x$observations,o2);x$links <- rbind(x$links,data.frame(observation_id="o2",entity_id="e1"))
p2 <- x$peptide_measurements;p2$observation_id <- "o2";x$peptide_measurements <- rbind(x$peptide_measurements,p2)
r <- pm_sample_ratios(x)$ratios
stopifnot(nrow(r)==4L,length(unique(r$observation_id))==2L,length(unique(r$modifications))==2L)

# Reject ambiguous canonical mappings instead of silently choosing a first row.
for(name in c("observations","entities","links","samples","peptide_measurements","protein_measurements")) {
  x <- fixture();x[[name]] <- rbind(x[[name]],x[[name]][1,,drop=FALSE]);fails(pm_sample_ratios(x))
}
for(column in c("peptide_column","protein_column")) {
  x <- fixture();x$samples[[column]][2] <- x$samples[[column]][1];fails(pm_sample_ratios(x))
}
x <- fixture();x$links$observation_id <- "unknown";fails(pm_sample_ratios(x))
x <- fixture();x$peptide_measurements$observation_id[1] <- "unknown";fails(pm_sample_ratios(x))
x <- fixture();x$protein_measurements$entity_id[1] <- "unknown";fails(pm_sample_ratios(x))
x <- fixture();x$protein_measurements$sample_id[1] <- "unknown";fails(pm_sample_ratios(x))
x <- fixture();x$peptide_measurements$source_column[1] <- "other column";fails(pm_sample_ratios(x))
x <- fixture();x$spec$peptide$scale <- NULL;fails(pm_sample_ratios(x))
x <- fixture();x$peptide_measurements$log2_value[1] <- 1;fails(pm_sample_ratios(x))
x <- fixture();x$observations$excluded <- NA;fails(pm_sample_ratios(x))
x <- fixture();x$peptide_measurements$input_value[1] <- 0;fails(pm_sample_ratios(x))

# No observations and no entities remain a valid, typed empty result.
x <- fixture()
for(name in c("observations","entities","links","peptide_measurements","protein_measurements")) x[[name]] <- x[[name]][FALSE,,drop=FALSE]
out <- pm_sample_ratios(x)
stopifnot(nrow(out$ratios)==0L,out$summary$ratio_rows==0L,is.logical(out$ratios$usable_pair),is.numeric(out$ratios$normalized_ratio))
cat("Sample ratio tests passed.\n")
