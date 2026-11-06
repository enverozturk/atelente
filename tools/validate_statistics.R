# Reproducible statistical implementation checks and simulation stress scenarios.
# Usage: Rscript tools/validate_statistics.R OUTPUT_DIRECTORY [REPLICATES]
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1L || length(args)>2L) stop("Usage: Rscript tools/validate_statistics.R OUTPUT_DIRECTORY [REPLICATES]")
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
for(module in c('mapping','input_adapter','sample_ratios','statistics')) source(file.path(root,'R',paste0(module,'.R')),encoding='UTF-8')
outdir <- args[1]
replicates <- if(length(args)==2L) suppressWarnings(as.numeric(args[2])) else 40
if(length(replicates)!=1L || !is.finite(replicates) || replicates<1 || replicates!=floor(replicates)) stop("REPLICATES must be a positive integer.")
if(file.exists(outdir) && !dir.exists(outdir)) stop("OUTPUT_DIRECTORY must be a directory.")
if(dir.exists(outdir) && length(list.files(outdir,all.files=TRUE,no..=TRUE))) stop("Choose an absent or empty output directory.")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
if(!dir.exists(outdir)) stop("Could not create output directory.")
ns <- 8L;nf <- 300L;ng <- 30L
sample_ids <- paste0('S',seq_len(ns));condition <- rep(c('Control','Treatment'),each=4)
samples <- data.frame(Sample=sample_ids,Condition=condition,Use='Yes','Experimental unit'=paste0('U',seq_len(ns)),check.names=FALSE)
comparisons <- data.frame(Name='Treatment-Control',Treatment='Treatment',Control='Control')
quant_samples <- data.frame(Sample=sample_ids,Condition=condition,Use='Yes','Peptide column'=sample_ids,'Protein column'=sample_ids,check.names=FALSE)
spec <- list(observation_level='peptide_observation',link_separator=NULL,
  peptide=list(id='ID',sequence='Sequence',protein_links='Group',scale='log2',quantity_type='Synthetic measured peptide log2 intensity'),
  protein=list(id='ID',accessions='Members',entity_type='protein_group',scale='log2',quantity_type='Synthetic measured group log2 intensity'))
near <- function(actual,expected,label) {
  if(!isTRUE(all.equal(as.numeric(actual),as.numeric(expected),tolerance=1e-10,check.attributes=FALSE))) stop(label)
}
generate <- function(seed,scenario) {
  set.seed(seed)
  group <- rep(seq_len(ng),each=nf/ng)
  beta <- numeric(nf)
  if(startsWith(scenario,'mixed')) beta[seq_len(60)] <- rep(c(1.5,-1.5),each=30)
  baseline <- rnorm(ng,24,1.2)
  biological_noise <- matrix(rnorm(ng*ns,0,.4),ng,ns)
  protein_shift <- outer(as.numeric(seq_len(ng) %% 2L==0L)*1.5,as.numeric(condition=='Treatment'))
  true_protein <- baseline+biological_noise+protein_shift
  shared_error <- matrix(rnorm(ng*ns,0,if(scenario=='null_independent') 0 else .35),ng,ns)
  protein <- true_protein+shared_error
  sigma <- sqrt(4/rchisq(nf,4))*.45
  independent_error <- matrix(rnorm(nf*ns),nf,ns)*sigma
  # Measured peptide/protein errors are correlated; common biological protein
  # changes cancel. Remaining group denominator noise induces row dependence.
  peptide <- true_protein[group,]+rnorm(nf,-4,.8)+outer(beta,as.numeric(condition=='Treatment'))+
    .7*shared_error[group,]+independent_error
  response_complete <- peptide-protein[group,]
  missing <- matrix(FALSE,nf,ns)
  if(grepl('mcar',scenario)) missing <- matrix(runif(nf*ns)<.15,nf,ns)
  if(grepl('mnar',scenario)) {
    residual_error <- independent_error-.3*shared_error[group,]
    missing[,condition=='Treatment'] <- residual_error[,condition=='Treatment']<0
  }
  peptide[missing] <- NA_real_
  pep <- data.frame(ID=sprintf('obs_%04d',seq_len(nf)),Sequence='PEPTIDE',Group=sprintf('G_%03d',group),check.names=FALSE)
  prot <- data.frame(ID=sprintf('G_%03d',seq_len(ng)),Members=paste0('P',seq_len(ng),';P',seq_len(ng),'-2'),check.names=FALSE)
  for(j in seq_len(ns)) {pep[[sample_ids[j]]] <- peptide[,j];prot[[sample_ids[j]]] <- protein[,j]}
  data <- pm_import_quant_tables(pep,prot,spec,quant_samples)
  ratios <- pm_sample_ratios(data)$ratios
  list(ratios=ratios,beta=beta,response=response_complete,missing=missing,group=group)
}
audit <- function(case,result,reference=FALSE) {
  x <- case$response;x[case$missing] <- NA_real_
  nc <- rowSums(is.finite(x[,1:4]));nt <- rowSums(is.finite(x[,5:8]))
  means_c <- rowMeans(x[,1:4],na.rm=TRUE);means_t <- rowMeans(x[,5:8],na.rm=TRUE)
  means_c[!is.finite(means_c)] <- NA_real_;means_t[!is.finite(means_t)] <- NA_real_
  sse <- rowSums((x[,1:4]-means_c)^2,na.rm=TRUE)+rowSums((x[,5:8]-means_t)^2,na.rm=TRUE)
  df <- nc+nt-2
  expected_eligible <- nc>=2 & nt>=2 & is.finite(sse) & sse>0
  tab <- result$results
  stopifnot(nrow(tab)==nf,identical(tab$observation_id,sprintf('obs_%04d',seq_len(nf))),
    identical(tab$tested,expected_eligible),identical(tab$eligible_for_fit,expected_eligible),
    all(tab$n_control==nc),all(tab$n_treatment==nt),nrow(result$measurements)==nf*ns,
    sum(result$measurements$used_for_fit)==sum((nc+nt)[expected_eligible]))
  near(tab$mean_control_log2_ratio,means_c,'Control means')
  near(tab$mean_treatment_log2_ratio,means_t,'Treatment means')
  near(tab$log2_effect,means_t-means_c,'Known OLS coefficients')
  k <- which(expected_eligible)
  near(tab$residual_variance[k],sse[k]/df[k],'Residual variances')
  near(tab$residual_df[k],df[k],'Residual df')
  near(tab$standard_error_unscaled[k],sqrt(1/nc[k]+1/nt[k]),'SE factors')
  hand_t <- (means_t[k]-means_c[k])/sqrt(tab$posterior_variance[k]*(1/nc[k]+1/nt[k]))
  near(tab$moderated_t[k],hand_t,'Moderated t reconstruction')
  near(tab$p_value[k],2*pt(-abs(hand_t),df=tab$total_df[k]),'Student t tail probabilities')
  near(tab$adjusted_p_value[k],p.adjust(tab$p_value[k],method='BH'),'BH family reconstruction')
  stopifnot(all(tab$eb_pool_size==length(k)),all(tab$n_family_tests==length(k)),
    all(is.na(tab$p_value[!expected_eligible])),all(is.na(tab$adjusted_p_value[!expected_eligible])))
  if(reference) {
    d <- cbind(Control=as.numeric(condition=='Control'),Treatment=as.numeric(condition=='Treatment'))
    f <- limma::eBayes(limma::contrasts.fit(limma::lmFit(x[k,,drop=FALSE],d),c(-1,1)),trend=FALSE,robust=FALSE)
    near(tab$p_value[k],f$p.value[,1],'Direct limma reference')
    near(tab$posterior_variance[k],f$s2.post,'Direct limma posterior variance')
    near(tab$total_df[k],f$df.total,'Direct limma total df')
  }
  invisible(TRUE)
}
records <- list();count <- 0L
scenarios <- c('null_independent','null_shared','mixed_shared','null_mcar','null_mnar','mixed_mcar')
for(scenario in scenarios) {
  for(iteration in seq_len(replicates)) {
    seed <- 20261030L+iteration
    case <- generate(seed,scenario)
    before <- serialize(case$ratios,NULL)
    result <- pm_test_ratios(case$ratios,samples,comparisons,design='independent_units')
    stopifnot(identical(before,serialize(case$ratios,NULL)))
    audit(case,result,reference=iteration==1)
    tab <- result$results;tested <- tab$tested;is_null <- case$beta==0
    reject <- tested & tab$adjusted_p_value<.05;reject[is.na(reject)] <- FALSE
    nominal <- tested & tab$p_value<.05;nominal[is.na(nominal)] <- FALSE
    count <- count+1L
    records[[count]] <- data.frame(scenario=scenario,iteration=iteration,seed=seed,
      total=nf,tested=sum(tested),null_tested=sum(tested & is_null),alternative_tested=sum(tested & !is_null),
      nominal_null_rejections=sum(nominal & is_null),bh_rejections=sum(reject),bh_false_rejections=sum(reject & is_null),
      bh_true_rejections=sum(reject & !is_null),fdp=if(sum(reject)) sum(reject & is_null)/sum(reject) else 0,
      any_false_rejection=any(reject & is_null),null_effect_bias=mean(tab$log2_effect[tested & is_null]),
      missing_cells=sum(case$missing),used_measurements=sum(result$measurements$used_for_fit))
    if(iteration %% 10==0 || iteration==replicates) cat(scenario,iteration,'/',replicates,'\n')
  }
}
runs <- do.call(rbind,records)
write.csv(runs,file.path(outdir,'simulation_runs.csv'),row.names=FALSE,na='')
summary <- lapply(split(runs,runs$scenario),function(x) {
  any_ci <- binom.test(sum(x$any_false_rejection),nrow(x))$conf.int
  rates <- x$nominal_null_rejections/x$null_tested
  rate_ci <- if(nrow(x)>1L) mean(rates)+c(-1,1)*qt(.975,nrow(x)-1)*sd(rates)/sqrt(nrow(x)) else c(NA_real_,NA_real_)
  list(scenario=x$scenario[1],replicates=nrow(x),features_per_replicate=nf,units_per_condition=4,
    mean_tested=mean(x$tested),mean_nominal_null_rejection_rate=mean(rates),
    nominal_rate_mc95_interval=unname(rate_ci),mean_false_discovery_proportion=mean(x$fdp),
    probability_any_false_bh_rejection=mean(x$any_false_rejection),any_false_bh_binomial95_interval=unname(any_ci),
    mean_true_discoveries=mean(x$bh_true_rejections),mean_null_effect_bias=mean(x$null_effect_bias),
    mean_missing_cells=mean(x$missing_cells))
})
# Broad pipeline sanity gates, fixed before inspecting outcomes. Calibration is
# reported with uncertainty rather than promoted to universal FDR guarantees.
stopifnot(all(runs$tested[runs$scenario %in% c('null_independent','null_shared','mixed_shared')]==nf),
  mean(runs$bh_true_rejections[runs$scenario=='mixed_shared'])>15,
  mean(runs$null_effect_bias[runs$scenario=='null_mnar'])>
    mean(runs$null_effect_bias[runs$scenario=='null_mcar'])+.1)
report <- list(result='PASS',meaning='Implementation checks passed; inspect scenario-specific calibration and the MNAR limitation separately.',seed_base=20261030L,iterations_per_scenario=replicates,limma_version=as.character(packageVersion('limma')),simulations=nrow(runs),
  candidate_hypotheses=nrow(runs)*nf,independent_formula_checks='Means, residual variance/df, SE, moderated t tails, BH, no source mutation and measurement counts on every run; direct limma reference first run per scenario.',
  design='Four independent units per condition. Correlated peptide/protein errors and shared group denominator noise; half the groups have a protein-only 1.5 log2 shift that should cancel.',
  mnar_limit='Deliberate condition-dependent censoring of negative errors in Treatment is an adverse stress scenario. No-imputation does not repair this selection bias.',
  scope='Implementation and controlled simulation validation only; no universal FDR guarantee, paired/batch/heteroscedastic model validation or biological benchmark.',
  scenarios=summary)
jsonlite::write_json(report,file.path(outdir,'simulation_summary.json'),auto_unbox=TRUE,pretty=TRUE,digits=10,na='null')
cat(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE,digits=6),'\n')
