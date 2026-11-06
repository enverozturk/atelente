# Cache orchestration uses real providers/readers and a deterministic transport.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
for (f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE)) source(f,encoding="UTF-8")
td <- tempfile("reference_access_");dir.create(td)
cache <- file.path(td,"cache")
calls <- 0L
fields <- c("CLE0000001","A01.001","Test","test","Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","reference","P00001","4","Test","Test","0","NULL","1-20","MS","","synthetic")
quoted <- paste0("'",fields,"'");quoted[c(15,18,19)] <- fields[c(15,18,19)]
raw <- charToRaw(paste0(paste(quoted,collapse="\t"),"\t\r\n"))
fetch <- function(url,path) {
  calls <<- calls+1L
  if(grepl("ftp.ebi.ac.uk",url,fixed=TRUE)) {
    writeBin(raw,path);headers <- list("last-modified"="Wed, 22 Feb 2023 12:10:00 GMT")
  } else {
    entry <- list(primaryAccession="P00001",sequence=list(value="MACDEFGHIK"),features=list())
    jsonlite::write_json(list(results=list(entry)),path,auto_unbox=TRUE)
    headers <- list("x-uniprot-release"="TEST_01","x-uniprot-release-date"="1 January 2026","x-total-results"="1")
  }
  list(status_code=200L,headers=headers,url=url,retrieved_utc="2026-11-02T00:00:00Z")
}
offline <- function(...)stop("HTTP must not run")
fails <- function(expr,pattern) {
  e <- tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e),fixed=TRUE))
}
m <- pm_reference_resolve("merops",cache_dir=cache,fetch=fetch)
u <- pm_reference_resolve("uniprot",cache_dir=cache,accessions=c("P00001","Q00002"),fetch=fetch)
stopifnot(calls==2L,m$mode=="downloaded",u$mode=="downloaded")
mw <- pm_reference_resolve("merops",cache_dir=cache,fetch=offline)
uw <- pm_reference_resolve("uniprot",cache_dir=cache,accessions=c("Q00002","P00001"),fetch=offline)
us <- pm_reference_resolve("uniprot",cache_dir=cache,accessions="Q00002",fetch=offline)
stopifnot(mw$mode=="cache",uw$mode=="cache",identical(us$snapshot,u$snapshot))
stopifnot(pm_reference_info(m,"merops")$source_date=="Wed, 22 Feb 2023 12:10:00 GMT",
  pm_reference_info(u,"uniprot")$release=="TEST_01")
provided <- pm_reference_resolve("uniprot",list(snapshot=u$snapshot,provenance=u$provenance),
  cache_dir=file.path(td,"unused"),fetch=offline)
stopifnot(provided$mode=="supplied snapshot",!dir.exists(file.path(td,"unused")))
disabled <- pm_annotate_references(list(merops=list(),uniprot=list()),annotations=FALSE,
  cache_dir=file.path(td,"disabled"))
stopifnot(!disabled$reference_access$enabled,is.null(disabled$merops),is.null(disabled$uniprot),
  !dir.exists(file.path(td,"disabled")))
fails(pm_annotate_references(list(),FALSE,merops=list()),"conflicts")
fails(pm_reference_resolve("uniprot",cache_dir=cache,accessions=NA_character_,fetch=offline),"non-empty")
# Partial downloads and malformed responses never become reusable cache entries.
failed_cache <- file.path(td,"failed")
fails(pm_reference_resolve("merops",cache_dir=failed_cache,fetch=function(url,path){
  writeLines("partial",path);stop("simulated interruption")}),"simulated interruption")
stopifnot(!dir.exists(file.path(failed_cache,"merops","substrates-v1")))
recovered <- pm_reference_resolve("merops",cache_dir=failed_cache,fetch=fetch)
stopifnot(recovered$mode=="downloaded")
# The complete public workflow consumes the validated providers without changing
# the supplied reference sequence. Warm-cache runs must not call transport.
writeLines(c(">sp|P00001|TEST","MACDEFGHIK"),file.path(td,"input.fasta"))
write.csv(data.frame(Sequence="MACD","Master Protein Accessions"="P00001",check.names=FALSE),
  file.path(td,"input.csv"),row.names=FALSE)
cfg <- list(peptides="input.csv",fasta="input.fasta",sequence_column="Sequence",
  output_dir="first",excel=FALSE,cache_dir="cache",cleavage=list(enzyme="trypsin",lysine_blocked=FALSE))
saved_fetch <- pm_reference_fetch
pm_reference_fetch <- offline
manifest <- pm_run(cfg,config_base=td,template=file.path(root,"inst","templates","report_template.html"))
stopifnot(manifest$reference_access$enabled,
  manifest$reference_access$providers$uniprot$access=="cache",
  manifest$reference_access$providers$uniprot$release=="TEST_01")
report <- paste(readLines(file.path(td,"first/report.html"),warn=FALSE),collapse="\n")
stopifnot(grepl("Reference sources",report,fixed=TRUE),grepl("TEST_01",report,fixed=TRUE))
cfg$output_dir <- "explicit"
cfg$merops <- list(snapshot=m$snapshot,provenance=m$provenance)
cfg$uniprot <- list(snapshot=u$snapshot,provenance=u$provenance)
explicit <- pm_run(cfg,config_base=td,template=file.path(root,"inst","templates","report_template.html"))
stopifnot(explicit$reference_access$providers$merops$access=="supplied snapshot")
pm_reference_fetch <- saved_fetch
# Altered cached bytes fail closed, without fetching a silent replacement.
cat("tamper",file=m$snapshot,append=TRUE)
fails(pm_reference_resolve("merops",cache_dir=cache,fetch=offline),"integrity check failed")
# Even independently edited metadata cannot assert a new accession request.
meta_path <- file.path(dirname(u$snapshot),"cache_meta.json")
meta <- jsonlite::read_json(meta_path,simplifyVector=TRUE)
meta$requested_accessions <- c("P00001","Q00002","EXTRA")
jsonlite::write_json(meta,meta_path,auto_unbox=TRUE)
fails(pm_reference_resolve("uniprot",cache_dir=cache,accessions="EXTRA",fetch=offline),"coverage disagrees")
# Reference metadata is escaped and explicit opt-out remains visible.
info <- list(enabled=TRUE,providers=list(uniprot=pm_reference_info(u,"uniprot")))
info$providers$uniprot$release <- '<script>alert("x")</script>'
html <- pm_report_reference_html(info)
stopifnot(grepl("&lt;script&gt;",html,fixed=TRUE),!grepl("<script>",html,fixed=TRUE),
  grepl("annotations = FALSE",pm_report_reference_html(disabled$reference_access),fixed=TRUE),
  identical(pm_report_reference_html(NULL),""))
cat("Reference cache, offline reuse, integrity, failure and report provenance tests passed.\n")
