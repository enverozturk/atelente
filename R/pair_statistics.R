# Experimental pair-level inference, separate from the peptide/protein API.
# Rows are distinct directed measurement pairs; annotations are not replicates.
pm_test_pair_matrix <- function(response, samples, membership, control, treatment,
                                minimum=3L) {
  if(!is.matrix(response) || !is.numeric(response) || is.null(rownames(response)) ||
     is.null(colnames(response)) || anyDuplicated(rownames(response)) ||
     anyDuplicated(colnames(response)) || anyNA(rownames(response)) ||
     any(!nzchar(rownames(response))) || anyNA(colnames(response)) ||
     any(!nzchar(colnames(response)))) stop("Unique pair and sample matrix names are required.")
  if(any(is.infinite(response))) stop("Infinite ratios are not valid measurements.")
  required <- c("sample","condition","experimental_unit")
  if(!is.data.frame(samples) || !all(required %in% names(samples))) stop("Explicit sample metadata required.")
  for(k in required) if(!is.character(samples[[k]]) || anyNA(samples[[k]]) ||
     any(!nzchar(samples[[k]]))) stop("Sample metadata must be nonempty text.")
  if(anyDuplicated(samples$sample) || anyDuplicated(samples$experimental_unit))
    stop("Each sample must represent a distinct declared independent unit.")
  if(!setequal(samples$sample,colnames(response))) stop("Matrix and metadata sample identities differ.")
  if(length(control)!=1L || length(treatment)!=1L || is.na(control) || is.na(treatment) ||
     identical(control,treatment) || !all(c(control,treatment) %in% samples$condition))
    stop("Declare distinct observed control and treatment conditions.")
  if(length(minimum)!=1L || is.na(minimum) || minimum<3L || minimum!=as.integer(minimum))
    stop("This pair model requires an integer minimum >=3.")
  if(!is.data.frame(membership) || !all(c("measurement_pair_id","family") %in% names(membership)))
    stop("Explicit pair/family membership required.")
  for(k in c("measurement_pair_id","family")) if(!is.character(membership[[k]]) ||
     anyNA(membership[[k]]) || any(!nzchar(membership[[k]]))) stop("Invalid membership identifiers.")
  if(any(!membership$measurement_pair_id %in% rownames(response))) stop("Unknown membership pair.")
  membership <- unique(membership[,c("measurement_pair_id","family"),drop=FALSE])
  if(!nrow(membership)) stop("No family members supplied.")
  metadata <- samples[match(colnames(response),samples$sample),,drop=FALSE]
  keep <- metadata$condition %in% c(control,treatment)
  metadata <- metadata[keep,,drop=FALSE]; response <- response[,keep,drop=FALSE]
  csel <- metadata$condition==control; tsel <- metadata$condition==treatment
  if(sum(csel)<minimum || sum(tsel)<minimum) stop("Insufficient independent units for declared minimum.")
  if(!requireNamespace("limma",quietly=TRUE)) stop("limma is required.")
  output <- list(); warnings <- list()
  for(family in unique(membership$family)) {
    ids <- membership$measurement_pair_id[membership$family==family]
    y <- response[ids,,drop=FALSE]
    nc <- rowSums(is.finite(y[,csel,drop=FALSE])); nt <- rowSums(is.finite(y[,tsel,drop=FALSE]))
    mc <- rowMeans(y[,csel,drop=FALSE],na.rm=TRUE); mt <- rowMeans(y[,tsel,drop=FALSE],na.rm=TRUE)
    df <- nc+nt-2L
    rss <- rowSums((y[,csel,drop=FALSE]-mc)^2,na.rm=TRUE)+rowSums((y[,tsel,drop=FALSE]-mt)^2,na.rm=TRUE)
    variance <- rss/df; variance[df<=0] <- NA_real_
    status <- rep("eligible",length(ids))
    status[nc<minimum | nt<minimum] <- "insufficient_replicates"
    status[status=="eligible" & (!is.finite(variance) | !is.finite(mt-mc))] <- "numeric_failure"
    status[status=="eligible" & variance==0] <- "zero_residual_variance"
    eligible <- which(status=="eligible")
    z <- data.frame(measurement_pair_id=ids,family=family,comparison=paste0(treatment,"_vs_",control),
      n_control=nc,n_treatment=nt,mean_control=mc,mean_treatment=mt,
      log2_effect=mt-mc,residual_df=df,residual_variance=variance,
      eligible_for_fit=status=="eligible",status=status,stringsAsFactors=FALSE,row.names=NULL)
    for(k in c("standard_error_unscaled","moderated_t","p_value","adjusted_p_value",
               "posterior_variance","prior_variance","prior_df","total_df")) z[[k]] <- NA_real_
    z$prior_df_infinite <- NA; z$tested <- FALSE
    z$eb_pool_size <- 0L; z$n_family_tests <- 0L
    warns <- character()
    if(length(eligible)<2L) z$status[eligible] <- "insufficient_features_for_moderation" else {
      fit <- withCallingHandlers({
        # Closed-form OLS for disjoint group indicators. These are exactly the
        # lmFit + contrasts.fit sufficient statistics, including row-specific
        # missingness. Avoid thousands of redundant row-wise QR decompositions.
        f <- list(coefficients=matrix((mt-mc)[eligible],ncol=1L),
                  stdev.unscaled=matrix(sqrt(1/nc[eligible]+1/nt[eligible]),ncol=1L),
                  sigma=sqrt(variance[eligible]),df.residual=df[eligible])
        limma::eBayes(f,trend=FALSE,robust=FALSE)
      },warning=function(w) {warns <<- unique(c(warns,conditionMessage(w)));invokeRestart("muffleWarning")})
      z$log2_effect[eligible] <- fit$coefficients[,1]
      z$residual_df[eligible] <- fit$df.residual; z$residual_variance[eligible] <- fit$sigma^2
      z$standard_error_unscaled[eligible] <- fit$stdev.unscaled[,1]
      z$moderated_t[eligible] <- fit$t[,1]; z$p_value[eligible] <- fit$p.value[,1]
      z$posterior_variance[eligible] <- fit$s2.post
      z$prior_variance[eligible] <- rep(fit$s2.prior,length.out=length(eligible))
      prior <- rep(fit$df.prior,length.out=length(eligible))
      z$prior_df_infinite[eligible] <- is.infinite(prior) & prior>0
      prior[!is.finite(prior)] <- NA_real_; z$prior_df[eligible] <- prior
      z$total_df[eligible] <- fit$df.total; z$eb_pool_size <- length(eligible)
      valid <- is.finite(z$p_value[eligible]) & z$p_value[eligible]>=0 & z$p_value[eligible]<=1 &
        is.finite(z$moderated_t[eligible]) & is.finite(z$posterior_variance[eligible]) &
        z$posterior_variance[eligible]>0 & is.finite(z$total_df[eligible]) & z$total_df[eligible]>0
      z$tested[eligible] <- valid
      z$status[eligible] <- ifelse(valid,"tested","model_numeric_failure")
      z$p_value[eligible[!valid]] <- NA_real_
      tested <- eligible[valid]
      z$adjusted_p_value[tested] <- p.adjust(z$p_value[tested],"BH")
      z$n_family_tests <- length(tested)
    }
    for(k in names(z)[vapply(z,is.numeric,logical(1))]) z[[k]][!is.finite(z[[k]])] <- NA_real_
    output[[family]] <- z; warnings[[family]] <- warns
  }
  list(results=do.call(rbind,output),warnings=warnings,
       settings=list(minimum=minimum,method="OLS limma",trend=FALSE,robust=FALSE,
                     limma_version=as.character(utils::packageVersion("limma"))))
}

# Independent arithmetic checks; does not certify model calibration.
pm_check_pair_arithmetic <- function(fitted,response,samples,control,treatment,tolerance=1e-8) {
  r <- fitted$results
  y <- response[r$measurement_pair_id,,drop=FALSE]
  cond <- samples$condition[match(colnames(y),samples$sample)]
  a <- y[,cond==control,drop=FALSE]; b <- y[,cond==treatment,drop=FALSE]
  na <- rowSums(is.finite(a)); nb <- rowSums(is.finite(b))
  ma <- rowMeans(a,na.rm=TRUE); mb <- rowMeans(b,na.rm=TRUE)
  j <- which(r$tested)
  close <- function(x,y) stopifnot(isTRUE(all.equal(unname(x),unname(y),tolerance=tolerance,check.attributes=FALSE)))
  if(length(j)) {
    close(r$log2_effect[j],(mb-ma)[j])
    close(r$residual_df[j],(na+nb-2)[j])
    s2 <- (rowSums((a-ma)^2,na.rm=TRUE)+rowSums((b-mb)^2,na.rm=TRUE))/(na+nb-2)
    close(r$residual_variance[j],s2[j])
    close(r$standard_error_unscaled[j],sqrt(1/na[j]+1/nb[j]))
    close(r$moderated_t[j],r$log2_effect[j]/sqrt(r$posterior_variance[j])/r$standard_error_unscaled[j])
    close(r$p_value[j],2*pt(-abs(r$moderated_t[j]),df=r$total_df[j]))
  }
  for(f in unique(r$family)) {
    q <- which(r$family==f & r$tested)
    close(r$adjusted_p_value[q],p.adjust(r$p_value[q],"BH"))
    stopifnot(!anyDuplicated(r$measurement_pair_id[r$family==f]))
  }
  invisible(TRUE)
}
