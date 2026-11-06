script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
if(.Platform$OS.type=="windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8")))
for(module in c("mapping","report","input_adapter","sample_ratios","normalization_io","statistics","statistics_io","report_quantities"))
  source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
fails<-function(x) stopifnot(inherits(try(force(x),silent=TRUE),"try-error"))
td<-tempfile("pm_report_quantities_");dir.create(td)
peptides<-data.frame(ID=c("001","001E2","3/4","004","005","006"),Sequence=c("PEPTIDEMK","PEPTIDEMK","ACDEFG","KLMNPQ","ACDEF","PEPTID"),
  Modifications=c("","Oxidation [M8]","","","",""),Link=c("G1","G1","G1;G2","G1","G1","G1"),
  p1=c(0,1,0,3,2,1),p2=c(1,3,4,4,5,2),p3=c(2,4,6,3,7,4),p4=c(4,5,8,5,8,NA),check.names=FALSE)
proteins<-data.frame(ID=c("G1","G2"),Members=c("P1;P1-2","P2"),r1=c(1,2),r2=c(1,2),r3=c(1,2),r4=c(1,2))
source_path<-file.path(td,"peptides.csv");pm_write_csv(peptides,source_path);peptides<-pm_read_table(source_path)
fa<-file.path(td,"tiny.fasta");writeLines(c(">sp|P1|ONE","MPEPTIDEMKACDEFGKLMNPQACDEFPEPTIDR",">sp|P1-2|TWO","PEPTIDEMKACDEFG"),fa)
mapped<-pm_map(peptides,pm_read_fasta(fa),scope="all",sequence_column="Sequence")
samples<-data.frame(Sample=c("01","02","03","04"),Condition=c("Control","Control","Treatment","Treatment"),Use="Yes",
  "Peptide column"=paste0("p",1:4),"Protein column"=paste0("r",1:4),check.names=FALSE)
spec<-list(observation_level="peptide_observation",link_separator=";",
  peptide=list(id="ID",sequence="Sequence",modifications="Modifications",protein_links="Link",scale="log2",quantity_type="Synthetic"),
  protein=list(id="ID",accessions="Members",entity_type="protein_group",scale="linear",quantity_type="Synthetic"))
data<-pm_import_quant_tables(peptides,proteins,spec,samples);ratios<-pm_sample_ratios(data)
norm<-file.path(td,"ratios");pm_export_sample_ratios(data,ratios,norm,list(peptides=source_path),excel=FALSE)
stat_samples<-data.frame(Sample=samples$Sample,Condition=samples$Condition,Use="Yes","Experimental unit"=paste0("unit",1:4),check.names=FALSE)
comparisons<-data.frame(Name=c("T-C","C-T"),Treatment=c("Treatment","Control"),Control=c("Control","Treatment"))
stat<-pm_test_ratios(ratios$ratios,stat_samples,comparisons)
stats<-file.path(td,"stats");pm_export_statistics(stat,stats,list(ratios=file.path(norm,"sample_ratios.csv"),ratio_manifest=file.path(norm,"normalization_manifest.json")),excel=FALSE)
snapshot<-serialize(mapped,NULL)
before<-vapply(c(list.files(norm,full.names=TRUE),list.files(stats,full.names=TRUE)),function(p)digest::digest(file=p,algo="sha256"),character(1))
out<-pm_report_quantities(mapped,norm,stats)
stopifnot(identical(snapshot,serialize(mapped,NULL)),identical(out$observations$observation_id,c("001","001E2","3/4","004","005","006")),
  identical(out$observations$peptide_id,mapped$peptides$peptide_id),out$observations$sequence[1]==out$observations$sequence[2],
  is.na(out$observations$modifications[1]),out$observations$modifications[2]=="Oxidation [M8]",
  nrow(out$ratios)==28L,nrow(out$results)==14L,nrow(out$measurements)==56L,
  all(out$results$status[out$results$observation_id=="3/4"]=="ambiguous_baseline"),
  identical(out$ratios$peptide_input_value,ratios$ratios$peptide_input_value),
  all.equal(out$results$p_value,stat$results$p_value,tolerance=1e-12)==TRUE,
  identical(out$results$peptide_id, mapped$peptides$peptide_id[match(out$results$observation_id,peptides$ID)]))
after<-vapply(names(before),function(p)digest::digest(file=p,algo="sha256"),character(1));stopifnot(identical(before,after))
plain<-pm_report_quantities(mapped,norm)
stopifnot(nrow(plain$results)==0L,nrow(plain$measurements)==0L,nrow(plain$ratios)==28L,is.null(plain$settings$statistics))
# Mapping-specific exclusions remain visible without deleting source quantities.
excluded<-mapped;excluded$peptides$mapping_status[1]<-"excluded_contaminant"
e<-pm_report_quantities(excluded,norm);stopifnot(e$observations$mapping_status[1]=="excluded_contaminant",nrow(e$ratios)==28L)

# A normalization made from an existing mapping source export has two known
# prefix columns; their identities and every original field must still match.
prefixed<-data;prefixed$source_peptides<-data.frame(peptide_id=mapped$peptides$peptide_id,source_row=mapped$peptides$source_row,peptides,check.names=FALSE)
pnorm<-file.path(td,"prefixed");pm_export_sample_ratios(prefixed,ratios,pnorm,list(peptides=source_path),excel=FALSE)
stopifnot(nrow(pm_report_quantities(mapped,pnorm)$ratios)==28L)

# Even identical sequences with different modifications cannot be swapped;
# altered unrelated source fields also invalidate identity.
bad<-mapped;bad$source_peptides<-bad$source_peptides[c(2,1,3:6),];fails(pm_report_quantities(bad,norm))
bad<-mapped;bad$source_peptides$p1[1]<-"999";fails(pm_report_quantities(bad,norm))
bad<-mapped;bad$peptides$modifications[1]<-"Oxidation [M8]";fails(pm_report_quantities(bad,norm))
bad<-mapped;bad$peptides$modifications_available[1]<-FALSE;fails(pm_report_quantities(bad,norm))
bad<-mapped;bad$peptides$peptide_id[2]<-bad$peptides$peptide_id[1];fails(pm_report_quantities(bad,norm))
bad<-mapped;bad$peptides$source_row[2]<-bad$peptides$source_row[1];fails(pm_report_quantities(bad,norm))

copy_dir<-function(from,name) {dest<-file.path(td,name);dir.create(dest);stopifnot(all(file.copy(list.files(from,full.names=TRUE),dest)));dest}
change<-function(dir,manifest,name,fun,rehash=TRUE) {
  path<-file.path(dir,paste0(name,".csv"));x<-pm_read_table(path);pm_write_csv(fun(x),path)
  if(rehash) {m<-jsonlite::read_json(file.path(dir,manifest),simplifyVector=TRUE);m$csv_sha256[[name]]<-digest::digest(file=path,algo="sha256");
    jsonlite::write_json(m,file.path(dir,manifest),auto_unbox=TRUE,pretty=TRUE,na="null",null="null")}
}
tamper<-copy_dir(norm,"tamper");change(tamper,"normalization_manifest.json","sample_ratios",function(x){x$log2_ratio[1]<-"999";x},FALSE)
fails(pm_report_quantities(mapped,tamper))
duplicate<-copy_dir(norm,"duplicate");change(duplicate,"normalization_manifest.json","sample_ratios",function(x)rbind(x,x[1,]))
fails(pm_report_quantities(mapped,duplicate))
wrong_obs<-copy_dir(norm,"wrong_observation");change(wrong_obs,"normalization_manifest.json","sample_ratios",function(x){x$observation_id[1]<-"missing";x})
fails(pm_report_quantities(mapped,wrong_obs))
bad_prefix<-copy_dir(pnorm,"bad_prefix");change(bad_prefix,"normalization_manifest.json","source_peptides",function(x){x$peptide_id[1]<-"different";x})
fails(pm_report_quantities(mapped,bad_prefix))
wrong_source<-copy_dir(norm,"wrong_source");change(wrong_source,"normalization_manifest.json","source_peptides",function(x){x$p1[1]<-"999";x})
fails(pm_report_quantities(mapped,wrong_source))
bad_stats<-copy_dir(stats,"bad_stats");change(bad_stats,"statistics_manifest.json","results",function(x)rbind(x,x[1,]))
fails(pm_report_quantities(mapped,norm,bad_stats))
bad_measure<-copy_dir(stats,"bad_measure");change(bad_measure,"statistics_manifest.json","measurements",function(x){x$peptide_input_value[1]<-"999";x})
fails(pm_report_quantities(mapped,norm,bad_measure))
missing_measure<-copy_dir(stats,"missing_measure");change(missing_measure,"statistics_manifest.json","measurements",function(x)x[-1,])
fails(pm_report_quantities(mapped,norm,missing_measure))
bad_p<-copy_dir(stats,"bad_p");change(bad_p,"statistics_manifest.json","results",function(x){x$p_value[1]<-"Inf";x})
fails(pm_report_quantities(mapped,norm,bad_p))
bad_label<-copy_dir(stats,"bad_label");change(bad_label,"statistics_manifest.json","measurements",function(x){x$analysis_condition[1]<-"Treatment";x})
fails(pm_report_quantities(mapped,norm,bad_label))
bad_unit<-copy_dir(stats,"bad_unit");change(bad_unit,"statistics_manifest.json","measurements",function(x){x$experimental_unit[1]<-"wrong";x})
fails(pm_report_quantities(mapped,norm,bad_unit))
bad_status<-copy_dir(stats,"bad_status");change(bad_status,"statistics_manifest.json","results",function(x){x$status[1]<-"ambiguous_baseline";x})
fails(pm_report_quantities(mapped,norm,bad_status))
fails(pm_report_quantities(mapped,pnorm,stats)) # The ratio manifest itself differs.
cleavage <- list(sites=data.frame(site_id=c("pair","one","triple"),neo_n_count=1L,neo_c_count=c(1,0,1)),
  evidence=data.frame(site_id=c("pair","pair","triple","one","pair","pair"),peptide_id=c("a","b","a","x","x","c"),
    role=c("N","spanning","C","N","C","C"),is_candidate=c(TRUE,FALSE,TRUE,TRUE,FALSE,TRUE)))
res <- data.frame(peptide_id=rep(c("a","b","c","x"),2),comparison=rep(c("T1","T2"),each=4),
  tested=c(TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE),p_value=c(.01,.04,NA,.0001,.2,NA,NA,.0001))
paired <- pm_paired_bh(res,cleavage)
stopifnot(identical(paired$counts,list(T1=2L,T2=1L)),
  identical(paired$results$paired_adjusted_p_value,c(.02,.04,NA,NA,.2,NA,NA,NA)),
  identical(paired$results[names(res)],res))
cleavage$sites$neo_c_count[] <- 0
empty <- pm_paired_bh(res,cleavage)
stopifnot(all(is.na(empty$results$paired_adjusted_p_value)),all(unlist(empty$counts)==0))
cat("Report quantity identity, typed values, ambiguity and hash checks passed.\n")
