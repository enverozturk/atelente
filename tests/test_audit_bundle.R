script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
source(file.path(root,"R","audit_bundle.R"),encoding="UTF-8")
td<-tempfile("atelente_audit_");dir.create(td)
writeLines("id,value\n001,0\n002,",file.path(td,"measurements.csv"))
writeLines("legacy workbook fixture",file.path(td,"proteomaster_results.xlsx"))
writeLines("main workbook fixture",file.path(td,"atelente_results.xlsx"))
writeLines("report fixture",file.path(td,"report.html"))
dir.create(file.path(td,"quantities","nested"),recursive=TRUE)
dir.create(file.path(td,"statistics"))
writeLines("sample,ratio\n001,1.25",file.path(td,"quantities","nested","sample_ratios.csv"))
writeLines("id,p\n001,0.0123",file.path(td,"statistics","results.csv"))
hash<-function(p)digest::digest(file=p,algo="sha256")
expected<-hash(file.path(td,"measurements.csv"))
nested_expected<-hash(file.path(td,"quantities","nested","sample_ratios.csv"))
manifest<-list(csv_sha256=list(measurements=expected))
jsonlite::write_json(manifest,file.path(td,"run_manifest.json"),auto_unbox=TRUE)
original_manifest<-hash(file.path(td,"run_manifest.json"))
out<-pm_finalize_audit(td,manifest)
stopifnot(out$artifact_layout=="atelente_bundle_v1",out$audit_file=="atelente_audit.zip",
 file.exists(file.path(td,"atelente_results.xlsx")),file.exists(file.path(td,"report.html")),
 !file.exists(file.path(td,"measurements.csv")),
 hash(file.path(td,"audit","measurements.csv"))==expected,
 hash(file.path(td,"audit","run_manifest.json"))==original_manifest,
 hash(file.path(td,"run_manifest.json"))==original_manifest,
 hash(file.path(td,"atelente_audit.zip"))==out$audit_sha256)
listing<-utils::unzip(file.path(td,"atelente_audit.zip"),list=TRUE)$Name
stopifnot("audit/measurements.csv"%in%listing,"audit/README.txt"%in%listing,
 "audit/quantities/nested/sample_ratios.csv"%in%listing,"audit/statistics/results.csv"%in%listing,
 !dir.exists(file.path(td,"quantities")),!dir.exists(file.path(td,"statistics")),
 hash(file.path(td,"audit","quantities","nested","sample_ratios.csv"))==nested_expected,
 !any(grepl("atelente_results.xlsx|report.html",listing)))
index<-jsonlite::read_json(file.path(td,"audit","audit_files.json"))
for(name in names(index$sha256))stopifnot(hash(file.path(td,"audit",name))==index$sha256[[name]])
# Refuse a linked generated directory before traversing or moving its contents.
outside<-tempfile("audit_external_");dir.create(outside);writeLines("private fixture",file.path(outside,"input.csv"))
unsafe<-tempfile("audit_link_");dir.create(unsafe)
linked<-suppressWarnings(if(.Platform$OS.type=="windows")Sys.junction(outside,file.path(unsafe,"quantities")) else file.symlink(outside,file.path(unsafe,"quantities")))
if(isTRUE(linked)) {
  error<-tryCatch({pm_finalize_audit(unsafe,list());NULL},error=identity)
  stopifnot(inherits(error,"error"),grepl("Unsafe",conditionMessage(error)),
    file.exists(file.path(outside,"input.csv")),!file.exists(file.path(unsafe,"atelente_audit.zip")))
}
cat("External audit packaging: original bytes/manifests, archive inventory and hashes preserved.\n")
