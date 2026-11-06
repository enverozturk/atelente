# Exercise the public attachment CLI before packaging, using its verified fixture.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
source(file.path(root,"tests/test_report_refresh.R"),encoding="UTF-8")
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE)) source(f,encoding="UTF-8")
invoke_repack <- function(input,out) {
  log <- suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),
    shQuote(c(file.path(root,"tools/refresh_output_bundle.R"),input,out)),stdout=TRUE,stderr=TRUE))
  status <- attr(log,"status");if(is.null(status))status<-0L
  list(status=status,log=paste(log,collapse="\n"))
}
files_hash <- function(folder) {
  files <- list.files(folder,recursive=TRUE,full.names=TRUE)
  setNames(vapply(files,hash,character(1)),substring(files,nchar(folder)+2L))
}
# The cleavage template stores a compact intensity summary instead of samples.
cleaved <- payload(file.path(td,"compact_integrated/report.html"))
cleaved$cleavage <- pm_build_cleavage(mapped,pm_cleavage_profile(list(enzyme="trypsin",lysine_blocked=FALSE)),hash(fpath))$cleavage
stopifnot(nrow(cleaved$cleavage$sites)>0L,nrow(cleaved$cleavage$evidence)>0L)
pm_write_report(cleaved,file.path(root,"inst","templates","report_template.html"),file.path(td,"compact_integrated/report.html"))
for(name in c("integrated","quantity_only","compact_integrated")) {
  current <- file.path(td,name);original_payload<-payload(file.path(current,"report.html"))
  original_hashes <- files_hash(current)
  if(name=="compact_integrated") stopifnot(nrow(original_payload$quantities$intensity_summary)>0L)
  for(pass in 1:2) {
    destination <- file.path(td,paste0(name,"_repack_",pass))
    result <- invoke_repack(file.path(current,"report.html"),destination)
    if(result$status!=0L) stop(result$log)
    refreshed <- payload(file.path(destination,"report.html"))
    stopifnot(identical(refreshed$quantities,original_payload$quantities),
      identical(refreshed$matches,original_payload$matches),
      identical(refreshed$cleavage,original_payload$cleavage),
      identical(refreshed$reference_access,original_payload$reference_access))
    for(link in names(original_payload$quantity_files)) {
      old <- original_payload$quantity_files[[link]]
      new <- refreshed$quantity_files[[link]]
      stopifnot(identical(new,if(startsWith(old,"audit/"))old else paste0("audit/",old)),
        file.exists(file.path(destination,new)),
        hash(file.path(td,name,old))==hash(file.path(destination,new)))
    }
    index <- jsonlite::read_json(file.path(destination,"audit/audit_files.json"))
    for(link in unlist(refreshed$quantity_files)) {
      relative <- substring(link,nchar("audit/")+1L)
      stopifnot(index$sha256[[relative]]==hash(file.path(destination,link)))
    }
    zipped <- unzip(file.path(destination,"atelente_audit.zip"),list=TRUE)$Name
    stopifnot(all(unlist(refreshed$quantity_files) %in% zipped))
    current <- destination
  }
  stopifnot(identical(original_hashes,files_hash(file.path(td,name))))
}
# A failure after copying must not leave an apparently completed output folder.
broken <- file.path(td,"broken");dir.create(broken)
broken_payload<-original_payload
broken_payload$quantity_files$ratios<-"quantities/missing.csv"
pm_write_report(broken_payload,file.path(root,"inst","templates","report_template.html"),file.path(broken,"report.html"))
failed <- file.path(td,"failed_repack")
before <- list.files(td,all.files=TRUE)
result <- invoke_repack(file.path(broken,"report.html"),failed)
stopifnot(result$status!=0L,!dir.exists(failed),identical(before,list.files(td,all.files=TRUE)))
cat("Quantity/statistics attachment repack: values, links, audit hashes, repeated refresh, source preservation and failure cleanup passed.\n")
