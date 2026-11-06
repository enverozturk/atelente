script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root,"R/uniprot.R"),encoding="UTF-8")
source(file.path(root,"R/reference_uniprot.R"),encoding="UTF-8")

# This standalone regression deliberately runs in the locale that previously
# truncated the UTF-8 gzip JSON fixture. No UTF-8 system locale is required.
stopifnot(nzchar(Sys.setlocale("LC_CTYPE","C")))
td <- tempfile("uniprot_download_");dir.create(td)
entry <- function(id,aliases=list()) list(primaryAccession=id,secondaryAccessions=aliases,
  entryAudit=list(entryVersion=4,sequenceVersion=2),sequence=list(value="MACDEFGHIK"),
  features=list(list(type="Signal",description="Literal \u03b1 <tag> & text",location=list(
    start=list(value=1,modifier="EXACT"),end=list(value=4,modifier="EXACT")),
    evidences=list(list(evidenceCode="ECO:0000269",source="PubMed",id="fixture")))))
json <- function(x) as.character(jsonlite::toJSON(x,auto_unbox=TRUE,null="null",na="null"))
write_raw <- function(s,path){con<-file(path,"wb");on.exit(close(con));writeBin(charToRaw(enc2utf8(s)),con)}
headers <- function(release="TEST_01",date="1 January 2026",link=NULL,total=NULL) {
  h<-list("x-uniprot-release"=release,"x-uniprot-release-date"=date)
  if(!is.null(link))h$link<-link
  if(!is.null(total))h[["x-total-results"]]<-as.character(total)
  h
}
page <- function(entries=list(),h=headers(),status=200L,raw=NULL,url=NULL) list(entries=entries,headers=h,status=status,raw=raw,url=url)
transport <- function(pages) {
  env <- new.env();env$n<-0L;env$urls<-character()
  env$fetch <- function(url,path){
    env$n<-env$n+1L;env$urls<-c(env$urls,url)
    if(env$n>length(pages))stop("Unexpected request")
    p<-pages[[env$n]];if(inherits(p,"error"))stop(p)
    write_raw(if(is.null(p$raw))json(list(results=p$entries)) else p$raw,path)
    list(headers=p$headers,status_code=p$status,url=if(is.null(p$url))url else p$url,retrieved_utc="2026-01-02T03:04:05Z")
  }
  env
}
outdir <- function(name) file.path(td,name)
next_url <- "https://rest.uniprot.org/uniprotkb/search?query=fixture&cursor=opaque&format=json&size=500"
e1 <- entry("P00001");e2 <- entry("Q00002",list("P99998"))
mock <- transport(list(page(list(e1),headers(link=paste0('<',next_url,'>; rel="next"'),total=2)),page(list(e2),headers(total=2))))
result <- pm_fetch_uniprot_reference(c("P00001","P00001-2","Q00002","P99998","P99999","P00001"),outdir("good"),mock$fetch)
loaded <- pm_read_uniprot_snapshot(result$snapshot,result$provenance)
audit <- read.csv(file.path(outdir("good"),"request_audit.csv"),stringsAsFactors=FALSE)
stopifnot(mock$n==2L,identical(mock$urls[2],next_url),length(loaded$entries)==2L,
  identical(loaded$entries[[1]],jsonlite::fromJSON(json(e1),simplifyVector=FALSE)),
  identical(loaded$entries[[2]],jsonlite::fromJSON(json(e2),simplifyVector=FALSE)),
  loaded$provenance$requested_accessions==5L,loaded$provenance$release=="TEST_01",
  identical(sort(unlist(loaded$provenance$requested_accession_ids)),sort(c("P00001","P00001-2","Q00002","P99998","P99999"))),
  audit$status[audit$requested_accession=="P99998"]=="secondary_alias",
  audit$status[audit$requested_accession=="P00001-2"]=="not_returned",
  audit$status[audit$requested_accession=="P99999"]=="not_returned",
  grepl("accession:P00001-2",utils::URLdecode(mock$urls[1]),fixed=TRUE),
  !grepl("MACDEFGHIK",mock$urls[1],fixed=TRUE),
  loaded$provenance$requests[[1]]$response_sha256==digest::digest(file=file.path(outdir("good"),"responses/batch_000_page_000.json"),algo="sha256"))
# An empty mapped set is valid and never accesses HTTP.
inactive <- list(primaryAccession="P00003",entryType="Inactive",
  inactiveReason=list(inactiveReasonType="DELETED",deletedReason="Removed by source"))
im <- transport(list(page(list(inactive))))
ir <- pm_fetch_uniprot_reference("P00003",outdir("inactive"),im$fetch)
il <- pm_read_uniprot_snapshot(ir$snapshot,ir$provenance)
ia <- read.csv(file.path(outdir("inactive"),"request_audit.csv"),stringsAsFactors=FALSE)
stopifnot(ia$status=="inactive",identical(il$entries[[1]],inactive),
  "P00003" %in% unlist(il$provenance$inactive_accessions))
fixture <- list(fasta=data.frame(record_id=1L,accession="P00003",sequence="MACDEFGHIK"),
  matches=data.frame(record_id=1L,peptide_id="p1",start=2L,end=5L))
annotated <- pm_annotate_uniprot(fixture,ir$snapshot,ir$provenance)
stopifnot(annotated$uniprot$records$status=="sequence_unavailable",
  nrow(annotated$uniprot$terminal_annotations)==0L)
empty <- pm_fetch_uniprot_reference(character(),outdir("empty"),function(...)stop("No HTTP allowed"))
el <- pm_read_uniprot_snapshot(empty$snapshot,empty$provenance)
stopifnot(length(el$entries)==0L,el$provenance$records==0L,is.null(el$provenance$release),el$provenance$retrieval_status=="empty_request")
# Bounded batches preserve literal IDs and combine only consistent releases.
ids <- sprintf("Q%05d",1:101)
many <- transport(list(page(lapply(ids[1:100],entry),headers(total=100)),page(list(entry(ids[101])),headers(total=1))))
mr <- pm_fetch_uniprot_reference(rev(ids),outdir("batches"),many$fetch)
stopifnot(many$n==2L,length(pm_read_uniprot_snapshot(mr$snapshot,mr$provenance)$entries)==101L,
  grepl("Q00101",many$urls[2],fixed=TRUE),!grepl("Q00101",many$urls[1],fixed=TRUE))
fails <- function(name,pages,ids="P00001") {
  tr <- transport(pages);dest<-outdir(name)
  err <- try(pm_fetch_uniprot_reference(ids,dest,tr$fetch),silent=TRUE)
  stopifnot(inherits(err,"try-error"),!file.exists(file.path(dest,"snapshot_provenance.json")))
  invisible(as.character(err))
}
stopifnot(grepl("mixed-release",fails("mixed",list(page(list(e1),headers(link=paste0('<',next_url,'>; rel="next"'))),page(list(e2),headers("TEST_02")))),fixed=TRUE))
fails("no_release",list(page(list(e1),h=list("content-type"="application/json"))))
fails("bad_json",list(page(raw="not json")))
fails("missing_results",list(page(raw='{"wrong":[]}')))
fails("missing_sequence",list(page(list(list(primaryAccession="P00001")))))
fails("bad_sequence",list(page(list(list(primaryAccession="P00001",sequence=list(value="M AC"))))))
fails("http_error",list(page(status=503L)))
fails("failed_request",list(simpleError("Network unavailable")))
fails("foreign_link",list(page(list(e1),headers(link='<https://example.com/data>; rel="next"'))))
fails("foreign_redirect",list(page(list(e1),url="https://example.com/data")))
fails("wrong_endpoint",list(page(list(e1),headers(link='<https://rest.uniprot.org/uniprotkb/stream?query=x>; rel="next"'))))
fails("multiple_next",list(page(list(e1),headers(link=paste0('<',next_url,'>; rel="next", <',next_url,'>; rel="next"')))))
fails("cyclic",list(page(list(e1),headers(link=paste0('<',next_url,'>; rel="next"'))),page(list(e2),headers(link=paste0('<',next_url,'>; rel="next"')))))
fails("truncated",list(page(list(e1),headers(total=2))))
fails("count_changed",list(page(list(e1),headers(link=paste0('<',next_url,'>; rel="next"'),total=2)),page(list(e2),headers(total=3))))
fails("duplicate_page",list(page(list(e1),headers(link=paste0('<',next_url,'>; rel="next"'))),page(list(e1))))
skipped <- pm_fetch_uniprot_reference(c("P00001 OR *","custom_\u03b1"),outdir("unsupported"),function(...)stop("No HTTP allowed"))
sp <- pm_read_uniprot_snapshot(skipped$snapshot,skipped$provenance)$provenance
stopifnot(sp$records==0L,sp$retrieval_status=="no_queryable_accessions",
  setequal(unlist(sp$unsupported_accession_ids),c("P00001 OR *","custom_\u03b1")))
con <- file(file.path(outdir("unsupported"),"request_audit.csv"),"rb")
audit_text <- rawToChar(readBin(con,"raw",n=10000));close(con);Encoding(audit_text)<-"UTF-8"
stopifnot(grepl("custom_\u03b1",audit_text,fixed=TRUE))
# Fresh staging guards against accidentally mixing a prior partial or complete run.
stopifnot(inherits(try(pm_fetch_uniprot_reference("P00001",outdir("good"),mock$fetch),silent=TRUE),"try-error"))
cat("UniProt R retrieval: pagination, batches, full records, UTF-8, release consistency, missing/alias audit and failure guards passed.\n")
