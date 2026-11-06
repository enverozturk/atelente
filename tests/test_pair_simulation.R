source('R/pair_simulation.R',encoding='UTF-8')
condition <- c(rep('WT',6),rep('KO',5))
mu <- seq(10,20,length.out=100); s <- rep(0.5,100); protein <- rep(paste0('g',1:20),each=5)
g <- pm_simulate_pair_forms(mu,s,protein,condition,'mixed','none',7)
stopifnot(identical(g,pm_simulate_pair_forms(mu,s,protein,condition,'mixed','none',7)),
          sum(g$effect!=0)==10,!any(g$mask))
# Same denominator is the same generated measurement, not fresh noise per pair.
pair1 <- g$observed[1,]-g$observed[5,];pair2 <- g$observed[2,]-g$observed[5,]
stopifnot(max(abs(pair1-pair2-(g$observed[1,]-g$observed[2,])))<1e-12)
pure <- g$means+g$error
stopifnot(max(abs(pair1-(pure[1,]-pure[5,])))<1e-12)
zero <- pm_simulate_pair_forms(mu,s,protein,condition,'null','none',8)
stopifnot(all(zero$effect==0))
intensity <- pm_simulate_pair_forms(mu,s,protein,condition,'null','intensity',8)
stopifnot(all(diff(intensity$loss_probability[order(intensity$latent)])<=0),
          identical(is.na(intensity$observed),intensity$mask))
conditional <- pm_simulate_pair_forms(mu,s,protein,condition,'null','condition',8)
stopifnot(all(conditional$loss_probability[,condition=='WT']==0.05),
          all(conditional$loss_probability[,condition=='KO'] %in% c(0.05,0.50)))
cat('PASS simulation: deterministic truth, shared denominator, shared-factor cancellation and missingness\n')
