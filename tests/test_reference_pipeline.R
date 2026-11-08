# Complete-workflow dispatch, failure propagation and offline CLI integration.
# Optional PACKAGE_LIBRARY also checks the installed public entrypoint/assets.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
api <- new.env(parent=globalenv())
for (path in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))
  sys.source(path,envir=api)
td <- tempfile("reference_pipeline_");dir.create(td)
td <- normalizePath(td,winslash="/",mustWork=TRUE)
fa <- file.path(td,"reference.fasta")
writeLines(c(">sp|TEST_A|TEST","MACDEFGHIK"),fa)
peptides <- data.frame(Sequence=c("MACD","EFGHIK","CDEFGH"),
  "Master Protein Accessions"="TEST_A",check.names=FALSE)
write.csv(peptides,file.path(td,"peptides.csv"),row.names=FALSE)
cfg <- list(peptides="peptides.csv",fasta="reference.fasta",output_dir="result",
  sequence_column="Sequence",excel=FALSE,cleavage=list(enzyme="trypsin",lysine_blocked=FALSE))
template <- file.path(root,"inst","templates","report_template.html")
fails <- function(expr,pattern) {
  e <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e),fixed=TRUE))
}

# Use real input/mapping functions, intercept only access and export boundaries.
probe <- new.env(parent=api)
probe$pm_run <- api$pm_run;environment(probe$pm_run) <- probe
calls <- list();exported <- NULL
probe$pm_annotate_references <- function(result,annotations=TRUE,merops=NULL,uniprot=NULL,cache_dir=NULL) {
  calls[[length(calls)+1L]] <<- list(annotations=annotations,merops=merops,uniprot=uniprot,cache_dir=cache_dir)
  result
}
probe$pm_export <- function(result,output_dir,template,inputs,...) {
  exported <<- list(result=result,output_dir=output_dir,template=template,inputs=inputs,options=list(...))
  list(probe=TRUE)
}
stopifnot(isTRUE(probe$pm_run(cfg,td,template)$probe),isTRUE(calls[[1]]$annotations))
baseline <- api$pm_build_cleavage(api$pm_map(api$pm_read_table(file.path(td,"peptides.csv")),
  api$pm_read_fasta(fa),sequence_column="Sequence"),cfg$cleavage,digest::digest(file=fa,algo="sha256"))
stopifnot(identical(exported$result,baseline),identical(exported$output_dir,file.path(td,"result")),
  identical(exported$inputs$fasta,fa))
override <- cfg;override$cache_dir <- "cache"
override$merops <- list(snapshot="refs/m.csv.gz",provenance="refs/m.json")
override$uniprot <- list(snapshot="refs/u.json.gz",provenance="refs/u.json")
probe$pm_run(override,td,template)
stopifnot(isTRUE(calls[[2]]$annotations),identical(calls[[2]]$cache_dir,file.path(td,"cache")),
  identical(calls[[2]]$merops$snapshot,file.path(td,"refs/m.csv.gz")),
  identical(calls[[2]]$uniprot$provenance,file.path(td,"refs/u.json")))
bad <- cfg;bad$annotations <- "FALSE";fails(probe$pm_run(bad,td,template),"annotations must be TRUE or FALSE")
bad <- cfg;bad$merops <- list(snapshot="only");fails(probe$pm_run(bad,td,template),"merops requires")
before <- length(calls);bad <- cfg;bad$neo_spanning_samples <- "missing.xlsx"
fails(probe$pm_run(bad,td,template),"does not exist")
stopifnot(length(calls)==before)
probe$pm_annotate_references <- function(...) stop("reference acquisition failed",call.=FALSE)
exported <- NULL
fails(probe$pm_run(cfg,td,template),"reference acquisition failed")
stopifnot(is.null(exported),!dir.exists(file.path(td,"result")))

# Explicit opt-out uses the actual reference stage and never requires a cache.
offline <- cfg;offline$annotations <- FALSE
api$pm_run(offline,td,template)
stopifnot(file.exists(file.path(td,"result/report.html")),!dir.exists(file.path(td,"cache")))
manifest <- jsonlite::read_json(file.path(td,"result/run_manifest.json"),simplifyVector=TRUE)
stopifnot(identical(manifest$reference_access$enabled,FALSE),
          identical(manifest$reference_access$status,"disabled"))
cli <- offline;cli$output_dir <- "cli-result"
jsonlite::write_json(cli,file.path(td,"config.json"),auto_unbox=TRUE,pretty=TRUE)
log <- system2(file.path(R.home("bin"),"Rscript"),c(shQuote(file.path(root,"scripts","run.R")),"--config",shQuote(file.path(td,"config.json"))),stdout=TRUE,stderr=TRUE)
if (!is.null(attr(log,"status"))) stop(paste(log,collapse="\n"))
stopifnot(file.exists(file.path(td,"cli-result/report.html")))
for (name in c("peptides.csv","matches.csv","cleavage_sites.csv","cleavage_evidence.csv")) {
  a<-file.path(td,"result","audit",name);b<-file.path(td,"cli-result","audit",name)
  stopifnot(file.exists(a),file.exists(b),
    identical(readBin(a,"raw",file.info(a)$size),readBin(b,"raw",file.info(b)$size)))
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args)) {
  stopifnot(length(args)==1L)
  ns <- loadNamespace("atelente",lib.loc=normalizePath(args[1],mustWork=TRUE))
  installed <- offline;installed$output_dir <- "installed-result"
  ns$pm_run(installed,td)
  stopifnot(file.exists(file.path(td,"installed-result/report.html")))
}
cat("Reference pipeline defaults, overrides, explicit opt-out, failure and CLI checks passed.\n")
