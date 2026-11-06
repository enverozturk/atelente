args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1]), winslash="/")))
source(file.path(root,"R","mapping.R")); source(file.path(root,"R","cleavage.R")); source(file.path(root,"R","report.R"))
profile <- pm_cleavage_profile(list(enzyme="trypsin",lysine_blocked=FALSE))
blocked <- pm_cleavage_profile(list(enzyme="trypsin",lysine_blocked=TRUE,pre_digest_n_label="TMTpro"))
unknown <- pm_cleavage_profile(list(enzyme="unknown"))
stopifnot(pm_n_label("2xTMTpro [K1; K4]; 1xTMTpro [N-Term]",TRUE,"TMTpro"),
  !pm_n_label("1xTMTpro [K1]",TRUE,"TMTpro"),
  !pm_n_label("1xAcetyl [N-Term]; 1xTMTpro [K1]",TRUE,"TMTpro"),
  !pm_n_label("1xOtherTMTpro [N-Term]",TRUE,"TMTpro"),
  !pm_n_label("",TRUE,"TMTpro"),is.na(pm_n_label(NA,TRUE,"TMTpro")),
  is.na(pm_n_label("1xTMTpro [N-Term]",FALSE,"TMTpro")),
  is.na(pm_n_label("1xTMTpro [N-Term]",TRUE,NULL)))
stopifnot(pm_bond_context("MACDEFGH",1)$context_display=="---M|ACDE",
  pm_bond_context("MACDEFGH",7)$context_display=="DEFG|H---",
  pm_bond_context("MACDEFGH",4)$context_display=="MACD|EFGH",
  pm_bond_context("MACDEFGH",4)$context_complete,
  !pm_bond_context("MACDEFGH",1)$context_complete)
stopifnot(inherits(try(pm_bond_context("MACD",0),silent=TRUE),"try-error"),
  inherits(try(pm_bond_context("MACD",4),silent=TRUE),"try-error"),
  inherits(try(pm_cleavage_profile(list(enzyme="trypsin")),silent=TRUE),"try-error"))
stopifnot(!pm_terminal("AKAA",2,"N",profile,FALSE)$is_candidate,
  pm_terminal("AKAA",2,"N",blocked,FALSE)$is_candidate,
  pm_terminal("ARAA",2,"N",blocked,TRUE)$is_candidate,
  !pm_terminal("ARAA",2,"C",blocked,TRUE)$is_candidate,
  pm_terminal("ARPA",2,"N",blocked,FALSE)$is_candidate,
  is.na(pm_terminal("AAAA",2,"N",unknown,NA)$is_candidate),
  !pm_terminal("MAAA",0,"N",blocked,TRUE)$is_candidate,
  !pm_terminal("MAAA",4,"C",blocked,TRUE)$is_candidate,
  pm_terminal("MAAA",1,"N",blocked,FALSE)$processing_status=="possible_initiator_met_removal",
  pm_terminal("MAAA",1,"N",blocked,FALSE)$is_candidate)

# Two duplicate FASTA records, one isoform, repeat placements and distinct modifications.
td <- tempfile("cleavage_");dir.create(td)
fa <- file.path(td,"ref.fasta")
writeLines(c(">sp|P1|TEST GN=TEST","MACDEFGHIKLMACDEFGHIKL",">sp|P1|DUP","MACDEFGHIKLMACDEFGHIKL",
  ">sp|P1-2|ISO","MACDEFGHIKQ"),fa)
peptides <- data.frame("Annotated Sequence"=c("MACD","EFGHIK","CDEFGH","EFGHIK","ACD","LMACD"),
  "Master Protein Accessions"=rep("P1; P1-2",6),
  "Modifications"=c("","1xTMTpro [N-Term]","","1xOxidation [M1]","",""),
  "Abundance: F1: Sample"=c(100,200,NA,300,400,0),check.names=FALSE)
x <- pm_map(peptides,pm_read_fasta(fa)); before <- x$measurements
out <- pm_build_cleavage(x,blocked,digest::digest(file=fa,algo="sha256"))
s <- out$cleavage$sites;e <- out$cleavage$evidence;t <- out$cleavage$termini
stopifnot(identical(before,out$measurements),nrow(t)==2*nrow(x$matches),all(s$bond_left>0),
  !anyDuplicated(s$site_id),length(unique(s$record_id))==3L)
site <- s[s$record_id==1 & s$bond_left==4,]
ev <- e[e$site_id==site$site_id,]
stopifnot(nrow(site)==1,site$context_display=="MACD|EFGH",site$neo_n_count==2,
  site$neo_c_count==2,site$spanning_count==1,
  all(ev$start[ev$role=="N"]==5),all(ev$end[ev$role=="C"]==4),
  all(ev$start[ev$role=="spanning"]<=4 & ev$end[ev$role=="spanning"]>=5),
  length(unique(s$site_id[s$bond_left==4]))==3)
# Enzyme-compatible unlabeled N observation remains noncandidate at a site seeded by a labelled observation.
writeLines(c(">sp|P2|TEST","MAAREFGHIK"),file.path(td,"label.fasta"))
labels <- data.frame("Annotated Sequence"=c("EFGHIK","EFGHIK","MAAR","AREFG"),
  "Master Protein Accessions"=rep("P2",4),"Modifications"=c("1xTMTpro [N-Term]","","",""),check.names=FALSE)
y <- pm_build_cleavage(pm_map(labels,pm_read_fasta(file.path(td,"label.fasta"))),blocked,digest::digest(file=file.path(td,"label.fasta"),algo="sha256"))
ys <- y$cleavage$sites[y$cleavage$sites$bond_left==4,]
ye <- y$cleavage$evidence[y$cleavage$evidence$site_id==ys$site_id,]
stopifnot(ys$neo_n_count==1,ys$neo_c_count==0,ys$spanning_count==1,
  nrow(ye)==4,!ye$is_candidate[ye$peptide_id=="pep_0000002"],
  ye$classification[ye$peptide_id=="pep_0000002"]=="digestion_compatible")
protein <- data.frame("Annotated Sequence"="MACDEFGHIKLMACDEFGHIKL","Master Protein Accessions"="P1",check.names=FALSE)
z <- pm_build_cleavage(pm_map(protein,pm_read_fasta(fa)),blocked,digest::digest(file=fa,algo="sha256"))
stopifnot(nrow(z$cleavage$sites)==0,nrow(z$cleavage$evidence)==0)
pm_export(out,file.path(td,"export"),file.path(root,"inst","templates","report_template.html"),list(fasta=fa))
ratios <- read.csv(file.path(td,"export","neo_spanning_ratios.csv"), fileEncoding="UTF-8-BOM")
at_bond <- ratios[ratios$site_id == site$site_id, ]
stopifnot(nrow(at_bond)==4, all(at_bond$spanning_peptide_id=="pep_0000003"),
  all(at_bond$spanning_status=="missing"), all(is.na(at_bond$ratio)))
# Identical sequences in separate input rows retain their own quantities. Repeated
# placements of one observation at a bond do not create additional measurements.
quant <- data.frame("Abundance: WT"=c(2,4,8,16,0,NA),
  "Abundance: control"=c(1,NA,0,1,-1,1), "Abundances (Grouped): pooled"=1,
  "Abundances (Normalized): WT"=rep(100,6),check.names=FALSE)
ev <- rbind(ev,ev,transform(ev[ev$role=="spanning",],peptide_id="pep_0000002"))
rat <- pm_neo_spanning_ratios(ev,pm_measurements(quant))
wt <- rat[rat$sample=="WT" & rat$spanning_peptide_id=="pep_0000003",]
control <- rat[rat$sample=="control" & rat$spanning_peptide_id=="pep_0000003",]
stopifnot(nrow(wt)==4, all(rat$neo_peptide_id!=rat$spanning_peptide_id))
stopifnot(isTRUE(all.equal(sort(wt$ratio),c(.25,.5,2))),
  isTRUE(all.equal(sort(wt$log2_ratio),c(-2,-1,1))),
  all(control$spanning_status=="zero"), all(is.na(control$log2_ratio)),
  any(control$neo_status=="missing"), any(control$neo_status=="negative"),
  nrow(pm_neo_spanning_ratios(ye,pm_measurements(quant),"grouped"))==0,
  all(pm_neo_spanning_ratios(ev,pm_measurements(quant),"normalized")$ratio==1),
  identical(names(rat),names(pm_neo_spanning_ratios(z$cleavage$evidence,pm_measurements(quant)))))
stopifnot(all(file.exists(file.path(td,"export",c("cleavage_sites.csv","cleavage_evidence.csv","peptide_termini.csv","report.html")))))
cat("Cleavage tests passed.\n")

# The same pair across two samples per condition, with a fifth excluded channel.
pair <- wt[1, ]
cmp <- pair[rep(1, 5), ]
cmp$source_column <- paste0("channel", 1:5)
cmp$log2_ratio <- c(1, 3, -1, -1, 99)
cmp$neo_status <- cmp$spanning_status <- "observed"
design <- list(samples=data.frame(Use=c(rep("Yes",4),"No"), Condition=c("A","A","B","B","A"),
  "Peptide column"=cmp$source_column,check.names=FALSE),
  comparisons=data.frame(Name=c("A_vs_B","B_vs_A"),Treatment=c("A","B"),Control=c("B","A")))
effect <- pm_neo_spanning_contrasts(cmp, design)
stopifnot(identical(effect$log2_fc,c(3,-3)), all(effect$status=="quantified_in_both"),
  all(effect$treatment_usable_pairs==2), all(effect$control_samples==2))
cmp$neo_status[3:4] <- "missing"
cmp$log2_ratio[3:4] <- NA_real_
effect <- pm_neo_spanning_contrasts(cmp,design)
stopifnot(identical(effect$status,c("neo_quantified_only_in_treatment","neo_quantified_only_in_control")),
  all(is.na(effect$fold_change)),identical(effect$treatment_neo,c(2L,0L)))
cmp$spanning_status[4] <- "missing"
stopifnot(all(pm_neo_spanning_contrasts(cmp,design)$status=="insufficient_coquantification"))
cmp$spanning_status[4] <- "observed"
cmp$neo_status[4] <- "zero"
stopifnot(all(pm_neo_spanning_contrasts(cmp,design)$status=="insufficient_coquantification"),
  identical(names(pm_neo_spanning_contrasts(cmp[FALSE,],design)),names(effect)))
cat("Descriptive neo/spanning contrast tests passed.\n")
