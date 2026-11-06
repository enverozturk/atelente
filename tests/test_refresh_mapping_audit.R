# Upgrade legacy report metadata only from a verified saved mapping audit.
script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-dirname(dirname(normalizePath(script,winslash="/")))
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))source(f,encoding="UTF-8")
d<-jsonlite::read_json(file.path(root,"tests/fixtures/export_bundle.json"),simplifyVector=TRUE)
mapping<-list(schema="atelente-column-mapping-1",quantities=data.frame(
 source_column=paste("Original",seq_len(nrow(d$pair_view$samples))),
 canonical_column=d$pair_view$samples[["Peptide column"]],input_scale="log2"))
td<-tempfile("atelente_refresh_audit_");dir.create(td)
sha<-function(p)digest::digest(file=p,algo="sha256")
write_json<-function(x,p)jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,null="null",na="null",digits=NA)
make_source<-function(name,layout="audit",embedded=FALSE){
 dir<-file.path(td,name);dir.create(dir)
 p<-d;if(embedded)p$column_mapping_audit<-mapping
 pm_write_report(p,file.path(root,"inst","templates","report_template.html"),file.path(dir,"report.html"))
 m<-list(tool="Atelente")
 if(!embedded){
   subdir<-if(layout=="audit")file.path(dir,"audit") else dir
   dir.create(subdir,showWarnings=FALSE)
   audit<-file.path(subdir,"column_mapping_audit.json");write_json(mapping,audit)
   if(layout=="audit"){
     index<-file.path(subdir,"audit_files.json")
     write_json(list(schema="atelente_audit_v1",sha256=list(column_mapping_audit.json=sha(audit))),index)
     m$audit_index_sha256<-sha(index)
   } else m$output_sha256<-list(column_mapping_audit.json=sha(audit))
 }
 write_json(m,file.path(dir,"run_manifest.json"));dir
}
invoke<-function(dir,label){
 out<-file.path(td,paste0(label,"_refreshed"))
 text<-suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),
   shQuote(c(file.path(root,"tools/refresh_output_bundle.R"),file.path(dir,"report.html"),out)),stdout=TRUE,stderr=TRUE))
 status<-attr(text,"status");if(is.null(status))status<-0L
 list(status=status,text=paste(text,collapse="\n"),out=out)
}
decode<-function(file){
 html<-paste(readLines(file,encoding="UTF-8",warn=FALSE),collapse="\n")
 marker<-'<script type="application/json" id="data">';at<-regexpr(marker,html,fixed=TRUE)[1]
 tail<-substring(html,at+nchar(marker));end<-regexpr('</script>',tail,fixed=TRUE)[1]
 pm_report_decode_json(substr(tail,1,end-1L))
}
for(layout in c("audit","root")){
 source<-make_source(paste0("legacy_",layout),layout)
 r<-invoke(source,paste0("legacy_",layout));stopifnot(r$status==0L)
 p<-decode(file.path(r$out,"report.html"));original<-decode(file.path(source,"report.html"))
 stopifnot(identical(p$column_mapping_audit$quantities,mapping$quantities),
   identical(p$pair_view$ratios,original$pair_view$ratios),identical(p$pair_view$statistics,original$pair_view$statistics))
 settings<-as.data.frame(readxl::read_xlsx(file.path(r$out,"atelente_results.xlsx"),sheet="Settings"))
 ss<-settings[settings$Section=="Samples",,drop=FALSE]
 stopifnot(identical(ss$`Source column`,mapping$quantities$source_column),
   all(ss$`Input scale`=="log2"),all(ss$`Export scale`=="linear"))
}
source<-make_source("current",embedded=TRUE);r<-invoke(source,"current");stopifnot(r$status==0L)
stopifnot(identical(decode(file.path(r$out,"report.html"))$column_mapping_audit,
 decode(file.path(source,"report.html"))$column_mapping_audit))

saved_d<-d;saved_mapping<-mapping
d$pair_view<-NULL;mapping$quantities<-mapping$quantities[FALSE,,drop=FALSE]
source<-make_source("mapping_only");r<-invoke(source,"mapping_only");stopifnot(r$status==0L)
stopifnot(length(decode(file.path(r$out,"report.html"))$column_mapping_audit$quantities)==0L)
d<-saved_d;mapping<-saved_mapping

source<-make_source("tampered");cat(" ",file=file.path(source,"audit","column_mapping_audit.json"),append=TRUE)
r<-invoke(source,"tampered");stopifnot(r$status!=0L,grepl("audit SHA256 mismatch",r$text,fixed=TRUE),!dir.exists(r$out))
source<-make_source("unverified",layout="root");write_json(list(tool="Atelente"),file.path(source,"run_manifest.json"))
r<-invoke(source,"unverified");stopifnot(r$status!=0L,grepl("no verified SHA256",r$text,fixed=TRUE),!dir.exists(r$out))
source<-make_source("wrong_columns",layout="root");bad<-mapping;bad$quantities$canonical_column[1]<-"Unexpected column"
path<-file.path(source,"column_mapping_audit.json");write_json(bad,path)
write_json(list(output_sha256=list(column_mapping_audit.json=sha(path))),file.path(source,"run_manifest.json"))
r<-invoke(source,"wrong_columns");stopifnot(r$status!=0L,grepl("canonical columns do not match",r$text,fixed=TRUE),!dir.exists(r$out))
cat("Legacy audit/root metadata restoration, source/scale labels, unchanged ratios/tests, current payload, and tamper/unverified/misaligned refusal passed.\n")
