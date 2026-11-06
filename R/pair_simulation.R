# Frozen milestone-1 sensitivity generator. Not a PKD missingness estimator.
pm_simulate_pair_forms <- function(mu,sd,protein,condition,truth,missingness,seed) {
  stopifnot(length(mu)==length(sd),length(mu)==length(protein),all(is.finite(mu)),
            all(is.finite(sd)),all(sd>0),all(condition %in% c('WT','KO')),
            truth %in% c('null','mixed'),missingness %in% c('none','random','condition','intensity'))
  set.seed(seed)
  n <- length(mu); k <- length(condition)
  beta <- numeric(n)
  if(truth=='mixed') {
    altered <- sample.int(n,floor(0.10*n))
    beta[altered] <- sample(c(-1,1),length(altered),replace=TRUE)
  }
  proteins <- unique(protein)
  shared <- matrix(rnorm(length(proteins)*k,sd=0.25),length(proteins),k)[match(protein,proteins),,drop=FALSE]
  loading <- matrix(rep(rnorm(k,sd=0.15),each=n),n,k)
  error <- matrix(rnorm(n*k),n,k)*sd
  means <- matrix(mu,n,k)+outer(beta,as.numeric(condition=='KO'))
  latent <- means+shared+loading+error
  loss <- matrix(0,n,k)
  if(missingness=='random') loss[,] <- 0.10
  if(missingness=='condition') {
    loss[,] <- 0.05
    loss[,condition=='KO'] <- ifelse(latent[,condition=='KO',drop=FALSE]>means[,condition=='KO',drop=FALSE],0.50,0.05)
  }
  if(missingness=='intensity') loss <- plogis((as.numeric(quantile(latent,0.10))-latent)/0.75)
  mask <- matrix(runif(n*k),n,k)<loss
  observed <- latent; observed[mask] <- NA_real_
  list(observed=observed,latent=latent,effect=beta,mask=mask,means=means,
       error=error,shared=shared,loading=loading,loss_probability=loss)
}

pm_simulation_metrics <- function(fit,pair_truth,truth_mode,missingness,iteration,seed,missing_fraction) {
  r <- fit$results; ans <- list()
  for(family in unique(r$family)) {
    z <- r[r$family==family,,drop=FALSE]
    actual <- pair_truth[z$measurement_pair_id]
    stopifnot(!anyNA(actual))
    null <- actual==0; discoveries <- z$tested & !is.na(z$adjusted_p_value) & z$adjusted_p_value<0.05
    fp <- sum(discoveries & null); tp <- sum(discoveries & !null); total <- sum(discoveries)
    nnull <- sum(z$tested & null); nalt <- sum(z$tested & !null)
    ans[[family]] <- data.frame(truth=truth_mode,missingness=missingness,iteration=iteration,seed=seed,
      family=family,graph_pairs=nrow(z),tested=sum(z$tested),tested_null=nnull,tested_alternative=nalt,
      graph_alternative=sum(!null),discoveries=total,false_discoveries=fp,true_discoveries=tp,
      fdp=if(total) fp/total else 0,
      nominal_null_fraction=if(nnull) sum(z$tested & null & !is.na(z$p_value) & z$p_value<0.05)/nnull else NA_real_,
      power_testable=if(nalt) tp/nalt else NA_real_,
      recovery_all_alternatives=if(sum(!null)) tp/sum(!null) else NA_real_,
      missing_form_fraction=missing_fraction,warning_count=length(fit$warnings[[family]]))
  }
  do.call(rbind,ans)
}

pm_summarize_simulations <- function(runs) {
  groups <- split(runs,interaction(runs$truth,runs$missingness,runs$family,drop=TRUE))
  do.call(rbind,lapply(groups,function(z) {
    out <- z[1,c('truth','missingness','family'),drop=FALSE]
    out$repetitions <- nrow(z); out$mean_tested <- mean(z$tested)
    out$zero_tested_runs <- sum(z$tested==0)
    for(k in c('fdp','nominal_null_fraction','power_testable','recovery_all_alternatives','missing_form_fraction')) {
      x <- z[[k]][is.finite(z[[k]])]
      m <- if(length(x)) mean(x) else NA_real_
      se <- if(length(x)>1) sd(x)/sqrt(length(x)) else NA_real_
      out[[paste0(k,'_mean')]] <- m; out[[paste0(k,'_mcse')]] <- se
      out[[paste0(k,'_lo95')]] <- max(0,m-1.96*se)
      out[[paste0(k,'_hi95')]] <- min(1,m+1.96*se)
    }
    out$null_any_discovery <- NA_real_; out$null_any_lo95 <- NA_real_;out$null_any_hi95 <- NA_real_
    if(z$truth[1]=='null') {
      b <- binom.test(sum(z$discoveries>0),nrow(z))
      out$null_any_discovery <- mean(z$discoveries>0)
      out$null_any_lo95 <- b$conf.int[1]; out$null_any_hi95 <- b$conf.int[2]
    }
    lo <- if(z$truth[1]=='null') out$null_any_lo95 else out$fdp_lo95
    out$clear_fdr_excess <- is.finite(lo) && lo>0.05
    out$interpretation <- if(out$clear_fdr_excess)
      if(z$missingness[1] %in% c('none','random')) 'baseline_calibration_failure' else 'stress_limitation' else
        'no_clear_excess_detected_not_proof_of_control'
    out
  }))
}
