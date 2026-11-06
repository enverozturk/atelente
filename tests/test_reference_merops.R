script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
source(file.path(root,"R/merops.R"),encoding="UTF-8")
source(file.path(root,"R/reference_merops.R"),encoding="UTF-8")
td <- tempfile("merops_fetch_");dir.create(td)
fields <- c("CLE0000001","A01.001","\tTest substrate's name","test",
  "Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","A reference\nwith tab\tand caf\u00e9 [one][two]",
  "P1","4","Test species","Test protease","0","NULL","1-20","MS","","synthetic")
encode <- function(x) {
  quoted <- paste0("'",x,"'")
  quoted[15] <- x[15];quoted[18] <- x[18];quoted[19] <- x[19]
  paste0(paste(quoted,collapse="\t"),"\t\r\n")
}
two <- fields;two[1] <- "CLE0000002";two[5] <- "Ac";two[13] <- "No references";two[15] <- "-1"
raw <- charToRaw(iconv(paste0(encode(fields),encode(two)),from="UTF-8",to="latin1"))
fetch <- function(url,path) {
  stopifnot(identical(url,pm_merops_source_url()))
  writeBin(raw,path)
  list(status_code=200L,headers=list("last-modified"="Wed, 22 Feb 2023 12:10:00 GMT"),
       url=url,retrieved_utc="2026-11-02T00:00:00Z")
}
out <- pm_fetch_merops_reference(file.path(td,"valid"),fetch)
loaded <- pm_read_merops_snapshot(out$snapshot,out$provenance)
r <- loaded$reference;p <- loaded$provenance
stopifnot(nrow(r)==2L,r$substrate_name[1]==fields[3],r$reference_raw[1]==fields[13],
  identical(r$p4_p4prime,c("ARNDCQEG","")),identical(r$canonical_octamer,c("true","false")),
  identical(r$reference_ids,c("one;two","")),r$source_column_18_raw[1]=="0",r$source_column_19_raw[1]=="",
  identical(r$p1_position,c("4","-1")),is.null(p$release_version),
  p$source_last_modified=="Wed, 22 Feb 2023 12:10:00 GMT",p$records==2L,
  p$source_sha256==digest::digest(raw,algo="sha256",serialize=FALSE))
fails <- function(expr) stopifnot(inherits(try(expr,silent=TRUE),"try-error"))
fails(pm_fetch_merops_reference(file.path(td,"valid"),fetch))
for (bad in list(charToRaw("<html>Not a substrate export</html>"),raw[-length(raw)][-1],
                charToRaw(paste0(encode(fields),encode(fields))),as.raw(c(1,0,2)))) {
  f <- function(url,path){writeBin(bad,path);list(status_code=200L,headers=list(),retrieved_utc="now")}
  stage <- tempfile(tmpdir=td)
  fails(pm_fetch_merops_reference(stage,f))
  stopifnot(!file.exists(file.path(stage,"snapshot_provenance.json")))
}
fails(pm_fetch_merops_reference(file.path(td,"transport"),function(url,path) stop("HTTP failure")))
fails(pm_fetch_merops_reference(file.path(td,"missing"),function(url,path) list(status_code=200L)))
fails(pm_fetch_merops_reference(file.path(td,"status"),function(url,path){writeBin(raw,path);list(status_code=503L)}))
cat("R-native MEROPS transport, parser, UTF-8 and provenance tests passed.\n")

# Optional local equivalence check: no network and no real files required by CI.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)) {
  stopifnot(length(args)==3L)
  local_fetch <- function(url,path) {
    stopifnot(file.copy(args[1],path))
    list(status_code=200L,headers=list(),url=url,retrieved_utc="local-equivalence-check")
  }
  real <- pm_fetch_merops_reference(file.path(td,"local_real"),local_fetch)
  a <- pm_read_merops_snapshot(real$snapshot,real$provenance)$reference
  b <- pm_read_merops_snapshot(args[2],args[3])$reference
  stopifnot(identical(names(a),names(b)),identical(a,b))
  cat("Full local MEROPS equivalence:",nrow(a),"records and",ncol(a),"columns identical.\n")
}
