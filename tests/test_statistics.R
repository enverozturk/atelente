args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1]),winslash="/")))
source(file.path(root,"R","sample_ratios.R"),encoding="UTF-8")
source(file.path(root,"R","statistics.R"),encoding="UTF-8")
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
near <- function(a,b,tolerance=1e-11) stopifnot(isTRUE(all.equal(a,b,tolerance=tolerance,check.attributes=FALSE)))
default_values <- rbind(c(-1,0,1,1,2,3),c(0,1,2,-1,0,1),c(-2,0,2,0,2,4),c(1,2,4,3,4,7))
fixture <- function(values=default_values,conditions=rep(c("Control","Treatment"),each=3)) {
  values <- as.matrix(values);stopifnot(ncol(values)==length(conditions))
  no <- nrow(values);ns <- ncol(values);oi <- rep(seq_len(no),each=ns);si <- rep(seq_len(ns),no)
  r <- pm_empty_sample_ratios()[rep(NA_integer_,no*ns),,drop=FALSE]
  r$observation_id <- paste0("o",oi);r$source_row <- oi+1L;r$sequence <- rep("PEPTIDE",length(oi))
  r$modifications <- ifelse(oi==2L,"Oxidation [M1]","");r$modifications_available <- TRUE
  r$entity_id <- paste0("e",oi);r$entity_type <- "protein_group";r$accessions <- paste0("P",oi,";P",oi,"-2")
  r$entity_source_row <- oi+1L;r$candidate_entity_count <- 1L;r$entity_status <- "matched"
  r$sample_id <- paste0("S",si);r$condition <- conditions[si]
  r$peptide_column <- paste0("p",si);r$protein_column <- paste0("r",si);r$peptide_scale <- "log2";r$protein_scale <- "log2"
  r$log2_ratio <- as.numeric(t(values));r$normalized_ratio <- 2^r$log2_ratio
  r$peptide_input_value <- r$peptide_log2_value <- r$log2_ratio;r$protein_input_value <- r$protein_log2_value <- 0
  r$peptide_original_value <- as.character(r$log2_ratio);r$protein_original_value <- "0"
  r$peptide_value_status <- ifelse(is.na(r$log2_ratio),"missing","observed");r$protein_value_status <- "observed"
  r$usable_pair <- !is.na(r$log2_ratio);r$ratio_status <- ifelse(r$usable_pair,"normalized","peptide_missing")
  r$ratio_representation <- ifelse(r$usable_pair,"linear_and_log2","unavailable")
  r$observation_excluded <- FALSE;r$entity_excluded <- FALSE;r$observation_exclusion_reason <- "";r$entity_exclusion_reason <- ""
  rownames(r) <- NULL
  samples <- data.frame(Sample=paste0("S",seq_len(ns)),Condition=conditions,Use="Yes",
    "Experimental unit"=paste0("independent ",seq_len(ns)),check.names=FALSE,stringsAsFactors=FALSE)
  comparisons <- data.frame(Name="Treatment vs Control",Treatment="Treatment",Control="Control",stringsAsFactors=FALSE)
  list(ratios=r,samples=samples,comparisons=comparisons)
}
run <- function(x) do.call(pm_test_ratios,x)

# Direct known group means, variance, sample size and model arithmetic.
x <- fixture();before <- serialize(x,NULL);out <- run(x);r <- out$results
stopifnot(identical(before,serialize(x,NULL)),all(r$status=="tested"),all(r$tested),all(r$n_control==3L),all(r$n_treatment==3L),
  all(r$eb_pool_size==4L),all(r$n_family_tests==4L),all(out$measurements$used_for_fit),
  nrow(out$measurements)==24L,nrow(r)==4L,length(unique(r$modifications))==2L,
  !any(grepl("significant|active_protease",names(r))))
near(r$mean_control_log2_ratio,rowMeans(default_values[,1:3]))
near(r$mean_treatment_log2_ratio,rowMeans(default_values[,4:6]))
near(r$log2_effect,c(2,-1,2,7/3));near(r$residual_df,rep(4,4));near(r$residual_variance,c(1,1,4,10/3))
near(r$standard_error_unscaled,rep(sqrt(2/3),4))
near(r$moderated_t,r$log2_effect/sqrt(r$posterior_variance*2/3))
near(r$p_value,2*stats::pt(-abs(r$moderated_t),df=r$total_df))
near(r$adjusted_p_value,stats::p.adjust(r$p_value,"BH"))
stopifnot(identical(out$samples,x$samples),identical(out$comparisons,x$comparisons),out$settings$limma_version==as.character(utils::packageVersion("limma")))

# Compare to a separately constructed limma fit and literal numeric contrast.
d <- cbind(C=c(1,1,1,0,0,0),T=c(0,0,0,1,1,1))
reference <- limma::eBayes(limma::contrasts.fit(limma::lmFit(default_values,d),matrix(c(-1,1),ncol=1)),trend=FALSE,robust=FALSE)
near(r$log2_effect,reference$coefficients[,1]);near(r$moderated_t,reference$t[,1]);near(r$p_value,reference$p.value[,1])
near(r$posterior_variance,reference$s2.post);near(r$total_df,reference$df.total)

# Multiple comparisons are separate families; reversing direction changes sign only.
x$comparisons <- rbind(x$comparisons,data.frame(Name="Control vs Treatment",Treatment="Control",Control="Treatment"))
r <- run(x)$results
near(r$log2_effect[1:4],-r$log2_effect[5:8]);near(r$p_value[1:4],r$p_value[5:8])
near(r$adjusted_p_value[1:4],r$adjusted_p_value[5:8]);stopifnot(all(r$n_family_tests==4L))

# Missing observations enter neither group means nor model residuals; no imputation.
v <- default_values;v[1,c(1,4)] <- NA_real_;x <- fixture(v);out <- run(x);r <- out$results
stopifnot(r$n_control[1]==2L,r$n_treatment[1]==2L,r$residual_df[1]==2L,r$status[1]=="tested",
  sum(out$measurements$used_for_fit)==22L)
near(r$log2_effect[1],2);near(r$residual_variance[1],.5);near(r$standard_error_unscaled[1],1)
v <- default_values;v[1,c(1,2)] <- NA_real_;out <- run(fixture(v))
stopifnot(out$results$status[1]=="insufficient_replicates",is.na(out$results$p_value[1]),
  out$results$n_control[1]==1L,is.finite(out$results$log2_effect[1]),all(out$results$eb_pool_size==3L),
  !any(out$measurements$used_for_fit[out$measurements$observation_id=="o1"]))

# A perfect shift with zero within-group variance stays descriptive.
v <- rbind(default_values,c(-4,-4,-4,4,4,4),rep(0,6),c(-1,0,1,1,2,3)*1e-10)
out <- run(fixture(v));r <- out$results
stopifnot(all(r$status[5:6]=="zero_residual_variance"),all(is.na(r$p_value[5:6])),r$log2_effect[5]==8,
  r$status[7]=="tested",r$residual_variance[7]>0,all(r$eb_pool_size==5L),all(r$n_family_tests==5L))

# No automatic unmoderated fallback when only one otherwise eligible feature remains.
out <- run(fixture(default_values[1,,drop=FALSE]));r <- out$results
stopifnot(r$status=="insufficient_features_for_moderation",r$eligible_for_fit,!r$tested,
  r$eligible_observation_count==1L,r$eb_pool_size==0L,r$n_family_tests==0L,is.na(r$p_value),!any(out$measurements$used_for_fit))
out <- run(fixture(matrix(c(0,0,0,2,2,2),nrow=1)))
stopifnot(out$results$status=="zero_residual_variance",!any(out$measurements$used_for_fit),all(is.na(out$results$prior_df_infinite)))

# Infinite prior df is valid and exported separately from missing numeric df.
# Whether limma estimates an infinite prior for two equal residual variances
# depends on the limma version (3.58 returns Inf; 3.62 returns a finite value),
# so the expected encoding is taken from limma's own squeezeVar on the same input.
out <- run(fixture(rbind(default_values[1,],default_values[2,])))
expected_prior <- limma::squeezeVar(out$results$residual_variance,out$results$residual_df)$df.prior
stopifnot(all(out$results$status=="tested"),all(is.finite(out$results$total_df)),all(out$results$total_df<=8))
if(all(is.infinite(expected_prior))) {
  stopifnot(all(out$results$prior_df_infinite),all(is.na(out$results$prior_df)),all(out$results$total_df==8))
} else {
  stopifnot(!any(out$results$prior_df_infinite),all(is.finite(out$results$prior_df)),
    isTRUE(all.equal(out$results$prior_df,rep(expected_prior,length.out=nrow(out$results)),tolerance=1e-10)))
}
stopifnot(!any(vapply(out$results,function(v) is.numeric(v)&&any(is.infinite(v)),logical(1))))

# Finite but extreme log ratios can overflow an effect or residual variance.
v <- rbind(default_values,rep(c(-1e308,1e308),each=3),c(-1e200,0,1e200,0,1e200,2e200))
x <- fixture(v)
overflow <- !is.finite(x$ratios$normalized_ratio);underflow <- x$ratios$normalized_ratio==0
x$ratios$ratio_representation[overflow] <- "log2_only_overflow"
x$ratios$ratio_representation[underflow] <- "log2_only_underflow"
x$ratios$normalized_ratio[overflow | underflow] <- NA_real_
out <- run(x)
stopifnot(all(out$results$status[5:6]=="numeric_range_exceeded"),all(is.na(out$results$p_value[5:6])),
  !any(vapply(out$results,function(v) is.numeric(v)&&any(is.infinite(v)),logical(1))),all(out$results$eb_pool_size==4L))

# Alternative baseline rows never join the empirical-Bayes pool or BH family.
x <- fixture();alt <- x$ratios[x$ratios$observation_id=="o1",];alt$entity_id <- "alternative";alt$accessions <- "P_other";alt$entity_source_row <- 99L
x$ratios$candidate_entity_count[x$ratios$observation_id=="o1"] <- alt$candidate_entity_count <- 2L
x$ratios <- rbind(x$ratios,alt);out <- run(x);r <- out$results
stopifnot(sum(r$status=="ambiguous_baseline")==2L,all(r$eb_pool_size==3L),all(r$n_family_tests==3L),
  !any(out$measurements$used_for_fit[out$measurements$observation_id=="o1"]))
unique_only <- fixture(default_values[2:4,,drop=FALSE]);ref <- run(unique_only)$results
near(r$p_value[r$tested],ref$p_value);near(r$posterior_variance[r$tested],ref$posterior_variance)

# Explicit exclusions and missing links stay visible without a test.
diagnostic <- function(x,observation,status) {
  pick <- x$ratios$observation_id==observation
  x$ratios$usable_pair[pick] <- FALSE;x$ratios$ratio_status[pick] <- status
  x$ratios$ratio_representation[pick] <- "unavailable";x$ratios$log2_ratio[pick] <- x$ratios$normalized_ratio[pick] <- NA_real_
  x
}
x <- diagnostic(fixture(),"o1","excluded_observation")
x$ratios$observation_excluded[x$ratios$observation_id=="o1"] <- TRUE
x$ratios$observation_exclusion_reason[x$ratios$observation_id=="o1"] <- "source_contaminant_flag"
stopifnot(run(x)$results$status[1]=="excluded_observation")
x <- diagnostic(fixture(),"o1","missing_protein_link");pick <- x$ratios$observation_id=="o1"
x$ratios$entity_id[pick] <- NA_character_;x$ratios$candidate_entity_count[pick] <- 0L
x$ratios$entity_status[pick] <- "missing_protein_link";x$ratios$entity_excluded[pick] <- NA
stopifnot(run(x)$results$status[1]=="missing_protein_link")
x <- diagnostic(fixture(),"o1","peptide_zero")
stopifnot(run(x)$results$status[1]=="insufficient_replicates")

# Source conditions can be filled when blank but cannot silently be overridden.
x <- fixture();x$ratios$condition <- "";stopifnot(all(run(x)$results$status=="tested"))
x <- fixture();x$samples$Condition[1] <- "Treatment";fails(run(x))
x <- fixture();x$ratios$condition[1] <- "other";fails(run(x))

# Unselected samples and third conditions are absent from a contrast's model.
v <- cbind(default_values,1000:1003,2000:2003)
x <- fixture(v,c(rep(c("Control","Treatment"),each=3),"Other","Other"))
# Avoid unrepresentable linear displays; canonical rows retain finite log2 ratios.
overflow <- !is.finite(x$ratios$normalized_ratio)
x$ratios$normalized_ratio[overflow] <- NA_real_;x$ratios$ratio_representation[overflow] <- "log2_only_overflow"
out <- run(x);near(out$results$p_value,run(fixture())$results$p_value)
stopifnot(!any(out$measurements$analysis_condition=="Other"))
x <- fixture();x$samples$Use[c(3,6)] <- "No";out <- run(x)
stopifnot(all(out$results$n_control==2L),all(out$results$n_treatment==2L),nrow(out$measurements)==16L)

# Exact metadata contracts reject accidental paired/batch/technical replication.
x <- fixture();x$samples[["Experimental unit"]][2] <- x$samples[["Experimental unit"]][1];fails(run(x))
x <- fixture();x$samples[["Experimental unit"]][4] <- x$samples[["Experimental unit"]][1];fails(run(x))
x <- fixture();x$samples$Batch <- "B1";fails(run(x))
x <- fixture();x$samples$Pair <- "donor";fails(run(x))
x <- fixture();fails(pm_test_ratios(x$ratios,x$samples,x$comparisons,design="paired"))
x <- fixture();x$samples$Use[c(2,3)] <- "No";fails(run(x))
x <- fixture();x$samples$Sample[2] <- x$samples$Sample[1];fails(run(x))
x <- fixture();x$samples$Sample[1] <- "unknown";fails(run(x))
x <- fixture();x$samples$Condition[1] <- "";fails(run(x))
x <- fixture();x$comparisons$Treatment <- "Control";fails(run(x))
x <- fixture();x$comparisons$Treatment <- "unknown";fails(run(x))
x <- fixture();x$comparisons <- rbind(x$comparisons,x$comparisons);fails(run(x))
x <- fixture();other <- x$comparisons;other$Name <- "repeated pair";x$comparisons <- rbind(x$comparisons,other);fails(run(x))

# Malformed ratio schemas cannot quietly change identity, candidate count or data.
x <- fixture();x$ratios <- rbind(x$ratios,x$ratios[1,]);fails(run(x))
x <- fixture();x$ratios$candidate_entity_count[1] <- 2L;fails(run(x))
x <- fixture();x$ratios$candidate_entity_count[1:6] <- 2L;fails(run(x))
x <- fixture();x$ratios$modifications[1] <- "changed identity";fails(run(x))
x <- fixture();x$ratios$accessions[1] <- "different entity";fails(run(x))
x <- fixture();x$ratios <- x$ratios[-1,];fails(run(x))
x <- fixture();x$ratios$log2_ratio[1] <- Inf;fails(run(x))
x <- fixture();x$ratios$usable_pair[1] <- FALSE;fails(run(x))
x <- fixture();x$ratios$ratio_status[1] <- "peptide_missing";fails(run(x))
x <- fixture();x$ratios$entity_excluded[1] <- TRUE;fails(run(x))
x <- fixture();x$ratios$normalized_ratio[1] <- 0;fails(run(x))
cat("Statistics tests passed.\n")
