script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
if(.Platform$OS.type=="windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8")))
for(module in c("mapping","report","input_adapter","sample_ratios","normalization_io","statistics","statistics_io"))
  source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
td<-tempfile("pm_report_refresh_");dir.create(td)
p<-data.frame(ID=c("001","001E2","003"),Sequence=c("PEPTIDEMK","PEPTIDEMK","ACDEFG"),Modifications=c("","Oxidation [M8]",""),Link="G1",
  p1=c(0,1,0),p2=c(1,3,4),p3=c(2,4,6),p4=c(4,5,8),check.names=FALSE)
p[["Leading razor protein"]]<-c("P1-2","P1","")
proteins<-data.frame(ID="G1",Members="P1;P1-2",r1=1,r2=1,r3=1,r4=1)
ppath<-file.path(td,"peptides.csv");pm_write_csv(p,ppath);p<-pm_read_table(ppath)
fpath<-file.path(td,"tiny.fasta");writeLines(c(">sp|P1|ONE","MPEPTIDEMKACDEFGR",">sp|P1-2|TWO","PEPTIDEMKACDEFG"),fpath)
mapped<-pm_map(p,pm_read_fasta(fpath),scope="all",sequence_column="Sequence")
mapped$reference_access <- list(enabled=FALSE,status="disabled",message="Annotations disabled by user (annotations = FALSE).")
mapping<-file.path(td,"mapping");pm_export(mapped,mapping,file.path(root,"inst","templates","report_template.html"),list(peptides=ppath,fasta=fpath),excel=FALSE)
samples<-data.frame(Sample=c("01","02","03","04"),Condition=c("Control","Control","Treatment","Treatment"),Use="Yes",
  "Peptide column"=paste0("p",1:4),"Protein column"=paste0("r",1:4),check.names=FALSE)
spec<-list(observation_level="peptide_observation",link_separator=NULL,
  peptide=list(id="ID",sequence="Sequence",modifications="Modifications",protein_links="Link",scale="log2",quantity_type="Synthetic"),
  protein=list(id="ID",accessions="Members",entity_type="protein_group",scale="linear",quantity_type="Synthetic"))
data<-pm_import_quant_tables(p,proteins,spec,samples);ratios<-pm_sample_ratios(data)
norm<-file.path(td,"ratios");pm_export_sample_ratios(data,ratios,norm,list(peptides=ppath),excel=FALSE,notes="Synthetic ratio fixture")
ss<-data.frame(Sample=samples$Sample,Condition=samples$Condition,Use="Yes","Experimental unit"=paste0("U",1:4),check.names=FALSE)
stat<-pm_test_ratios(ratios$ratios,ss,data.frame(Name="T-C",Treatment="Treatment",Control="Control"))
# Deliberately tiny export values exercise JSON presentation precision, not an
# additional scientific result. The existing statistics fixtures test fitting.
stat$results$p_value[1]<-1.23456789012345e-120;stat$results$adjusted_p_value[1]<-3.70370367037035e-120
stats<-file.path(td,"statistics");pm_export_statistics(stat,stats,
  list(ratios=file.path(norm,"sample_ratios.csv"),ratio_manifest=file.path(norm,"normalization_manifest.json")),excel=FALSE,notes="Synthetic export precision fixture")
hash<-function(p) digest::digest(file=p,algo="sha256")
all_inputs<-c(list.files(mapping,full.names=TRUE),list.files(norm,full.names=TRUE),list.files(stats,full.names=TRUE))
before<-vapply(all_inputs,hash,character(1))
config_path<-file.path(td,"report.json")
write_config<-function(output,source="mapping",statpath="statistics") {
  config<-list(mapping_dir=source,ratios_dir="ratios",output_dir=output)
  if(!is.null(statpath)) config$statistics_dir<-statpath
  jsonlite::write_json(config,config_path,auto_unbox=TRUE,pretty=TRUE)
}
run<-function() {
  executable<-file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
  log<-suppressWarnings(system2(executable,c(shQuote(file.path(root,"scripts","refresh_report.R")),"--config",shQuote(config_path)),stdout=TRUE,stderr=TRUE))
  status<-attr(log,"status");if(is.null(status)) status<-0L
  list(status=status,log=log)
}
payload<-function(path) {
  html<-paste(readLines(path,encoding="UTF-8",warn=FALSE),collapse="\n")
  marker<-'<script type="application/json" id="data">';start<-regexpr(marker,html,fixed=TRUE)[1]+nchar(marker)
  tail<-substring(html,start);end<-regexpr('</script>',tail,fixed=TRUE)[1]
  pm_report_decode_json(substr(tail,1,end-1L))
}
write_config("integrated");execution<-run();if(execution$status!=0L) stop(paste(execution$log,collapse="\n"))
out<-file.path(td,"integrated");embedded<-payload(file.path(out,"report.html"));original<-payload(file.path(mapping,"report.html"))
for(name in names(original)) stopifnot(identical(original[[name]],embedded[[name]]))
stopifnot(identical(embedded$peptides$leading_razor_accessions,c("P1-2","P1",NA_character_)))
stopifnot(nrow(embedded$quantities$results)==3L,nrow(embedded$quantities$ratios)==12L,
  identical(embedded$quantities$observations$observation_id,c("001","001E2","003")),
  abs(embedded$quantities$results$p_value[1]/stat$results$p_value[1]-1)<1e-12,
  abs(embedded$quantities$results$adjusted_p_value[1]/stat$results$adjusted_p_value[1]-1)<1e-12,
  embedded$quantities$settings$notes$quantity=="Synthetic ratio fixture",
  embedded$quantities$settings$notes$statistics=="Synthetic export precision fixture")
for(pair in list(c(mapping,out),c(norm,file.path(out,"quantities")),c(stats,file.path(out,"statistics")))) {
  inputs<-list.files(pair[1],pattern="\\.(csv|json|xlsx)$",full.names=TRUE)
  for(path in inputs) stopifnot(hash(path)==hash(file.path(pair[2],basename(path))))
}
stopifnot(identical(before,vapply(all_inputs,hash,character(1))))
manifest<-jsonlite::read_json(file.path(out,"report_manifest.json"),simplifyVector=TRUE)
stopifnot(identical(manifest$reference_access,mapped$reference_access))
stopifnot(manifest$stage=="integrated_quantity_report",manifest$mapping_report_sha256==hash(file.path(mapping,"report.html")),
  manifest$report_sha256==hash(file.path(out,"report.html")))
protected<-hash(file.path(out,"report.html"));stopifnot(run()$status!=0L,protected==hash(file.path(out,"report.html")))

# Quantity-only attachment has no fabricated comparison or statistical result.
write_config("quantity_only",statpath=NULL);stopifnot(run()$status==0L)
q<-payload(file.path(td,"quantity_only","report.html"))
stopifnot(length(q$quantities$results)==0L,is.null(q$quantities$settings$statistics),nrow(q$quantities$ratios)==12L)

# Reordered source rows must fail before creating an output directory, even
# when the same peptide sequence appears twice.
bad<-file.path(td,"reordered");dir.create(bad);stopifnot(all(file.copy(list.files(mapping,full.names=TRUE),bad)))
raw<-pm_read_table(file.path(bad,"source_peptides.csv"));pm_write_csv(raw[c(2,1,3),],file.path(bad,"source_peptides.csv"))
write_config("refused",source="reordered");stopifnot(run()$status!=0L,!dir.exists(file.path(td,"refused")))
# A previously integrated report is not recursively refreshed/copied.
write_config("refused_twice",source="integrated");stopifnot(run()$status!=0L,!dir.exists(file.path(td,"refused_twice")))
# Compact source layout: the same original bytes now live below audit/.
for(folder in c(mapping,norm,stats)) {
  audit<-file.path(folder,"audit");dir.create(audit)
  files<-list.files(folder,pattern="\\.(csv|json|xlsx)$",full.names=TRUE)
  stopifnot(all(file.rename(files,file.path(audit,basename(files)))))
}
write_config("compact_integrated");execution<-run();if(execution$status!=0L)stop(paste(execution$log,collapse="\n"))
compact<-payload(file.path(td,"compact_integrated","report.html"))
stopifnot(identical(embedded$quantities$results,compact$quantities$results),
  hash(file.path(mapping,"audit","source_peptides.csv"))==hash(file.path(td,"compact_integrated","audit","source_peptides.csv")))
cat("Report refresh CLI, source preservation, metadata, precision and refusal checks passed.\n")
