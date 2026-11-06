# The source archive must contain only the installable package (DESCRIPTION,
# NAMESPACE, R/, man/, inst/), even when the working repository holds local
# analysis inputs or results. Builds a copy of the package root that also
# contains synthetic canary folders and checks the archive members.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
td <- tempfile("atelente-build-exclusions-"); dir.create(td)
copy <- file.path(td,"atelente"); dir.create(copy)
members <- c("DESCRIPTION","NAMESPACE",".Rbuildignore","R","man","inst","tools","docs","tests")
members <- members[file.exists(file.path(root,members))]
stopifnot(all(file.copy(file.path(root,members),copy,recursive=TRUE)))
# Local folders a user might create while analysing data, plus an arbitrary name.
canaries <- c("data","work","outputs","atelente_runs","results","my_local_analysis")
for(folder in canaries) {
  dir.create(file.path(copy,folder))
  writeLines("EXCLUSION_CANARY",file.path(copy,folder,"EXCLUSION_CANARY.txt"))
}
writeLines("EXCLUSION_CANARY",file.path(copy,"local_note.txt"))
out <- file.path(td,"build"); dir.create(out)
old <- setwd(out); on.exit(setwd(old),add=TRUE)
log <- suppressWarnings(system2(file.path(R.home("bin"),"R"),c("CMD","build","--no-build-vignettes",shQuote(copy)),
  stdout=TRUE,stderr=TRUE))
archive <- list.files(out,pattern="^atelente_.*\\.tar\\.gz$",full.names=TRUE)
if(length(archive)!=1L) stop("R CMD build failed:\n",paste(log,collapse="\n"))
files <- sub("^atelente/","",utils::untar(archive,list=TRUE))
files <- files[nzchar(files) & !grepl("/$",files)]
demo <- paste0("inst/extdata/neutral_demo/", c("normalize.json", "peptides.csv", "proteins.csv", "samples.csv"))
allowed <- grepl("^(DESCRIPTION|NAMESPACE|R/[^/]+\\.R|man/[^/]+\\.Rd|inst/templates/[^/]+|inst/extdata/config/[^/]+\\.example\\.json|build/[^/]+)$",files) | files %in% demo
if(!all(allowed)) stop("Unexpected files in the source archive: ",paste(files[!allowed],collapse=", "))
if(any(grepl("CANARY|local_note",files,ignore.case=TRUE))) stop("Local canary files entered the archive.")
templates <- c("report_template.html","cleavage_report_template.html","pair_summary.js","report_payload.js","report_workbook.js")
configs <- paste0(c("config","maxquant","quantification","report","statistics"),".example.json")
stopifnot(all(paste0("inst/templates/",templates) %in% files),
  all(paste0("inst/extdata/config/",configs) %in% files),
  all(demo %in% files),
  identical(sort(sub("^R/","",files[startsWith(files,"R/")])),sort(list.files(file.path(root,"R"),pattern="\\.R$"))))
cat("Source archive contains the package, templates, configs and four synthetic demo files; local folders are excluded.\n")
