source('R/pair_statistics.R',encoding='UTF-8')
set.seed(91)
m <- data.frame(sample=paste0('s',1:11),condition=c(rep('WT',6),rep('KO',5)),experimental_unit=paste0('u',1:11))
y <- matrix(rnorm(60*11),60,11,dimnames=list(paste0('p',1:60),m$sample))
y[1:10,7:11] <- y[1:10,7:11]+2
y[11,1:4] <- NA; y[12,] <- 3
y[13,] <- c(rep(1,6),rep(3,5))
membership <- data.frame(measurement_pair_id=rownames(y),family=rep(c('A','B'),each=30))
original <- y
f <- pm_test_pair_matrix(y,m,membership,'WT','KO')
pm_check_pair_arithmetic(f,y,m,'WT','KO')
stopifnot(identical(y,original),f$results$status[11]=='insufficient_replicates',
          f$results$status[12]=='zero_residual_variance',f$results$status[13]=='zero_residual_variance')
# Duplicated annotations must not alter pooling or results.
d <- pm_test_pair_matrix(y,m,rbind(membership,membership[1:10,]),'WT','KO')
stopifnot(identical(f,d))
# Same measurement may annotate another family without multiplying that family.
multi <- rbind(membership,data.frame(measurement_pair_id=rownames(y)[1:12],family='C'))
fm <- pm_test_pair_matrix(y,m[nrow(m):1,],multi,'WT','KO')
pm_check_pair_arithmetic(fm,y,m,'WT','KO')
stopifnot(nrow(fm$results[fm$results$family=='C',])==12,
          isTRUE(all.equal(fm$results[fm$results$family!='C',],f$results)))
# Compare with a separate direct limma fit for the eligible first family.
idx <- which(f$results$family=='A' & f$results$tested)
design <- model.matrix(~factor(m$condition,levels=c('WT','KO')))
direct <- limma::eBayes(limma::lmFit(y[idx,],design),trend=FALSE,robust=FALSE)
stopifnot(isTRUE(all.equal(unname(direct$p.value[,2]),f$results$p_value[idx],tolerance=1e-10)))
reverse <- pm_test_pair_matrix(y,m,membership,'KO','WT')$results
stopifnot(isTRUE(all.equal(reverse$p_value,f$results$p_value,tolerance=1e-10)),
          isTRUE(all.equal(reverse$log2_effect,-f$results$log2_effect,tolerance=1e-10)))
# All missing and exactly 3/3 eligibility; unknown IDs and repeated animals fail.
y[14,] <- NA; y[15,c(4:6,10:11)] <- NA
z <- pm_test_pair_matrix(y,m,membership,'WT','KO')$results
stopifnot(z$status[14]=='insufficient_replicates',z$n_control[15]==3,z$n_treatment[15]==3,z$tested[15])
bad <- m;bad$experimental_unit[2] <- bad$experimental_unit[1]
stopifnot(inherits(try(pm_test_pair_matrix(y,bad,membership,'WT','KO'),silent=TRUE),'try-error'))
badm <- membership;badm$measurement_pair_id[1] <- 'unknown'
stopifnot(inherits(try(pm_test_pair_matrix(y,m,badm,'WT','KO'),silent=TRUE),'try-error'))
cat('PASS pair statistics: direct fit, arithmetic, BH, duplicates, missingness, reverse contrast and units\n')
