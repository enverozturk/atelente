script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-dirname(dirname(normalizePath(script,winslash="/")))
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))source(f)
td<-tempfile("identifier_scope_");dir.create(td)
calls<-character()
fetch<-function(url,path){
  calls<<-c(calls,url)
  jsonlite::write_json(list(results=list(list(primaryAccession="P00001",sequence=list(value="MACDEFGHIK"),features=list()))),path,auto_unbox=TRUE)
  list(status_code=200L,url=url,retrieved_utc="2026-11-05T00:00:00Z",headers=list(
    "x-uniprot-release"="TEST","x-uniprot-release-date"="test date","x-total-results"="1"))
}
ids<-c("P00001","CON__TEST","NP_000001.1","custom_name","TEST1")
r<-pm_reference_resolve("uniprot",cache_dir=file.path(td,"cache"),accessions=ids,fetch=fetch)
stopifnot(length(calls)==1L,grepl("accession:P00001",utils::URLdecode(calls),fixed=TRUE),
  !any(vapply(ids[-1],function(id)grepl(id,utils::URLdecode(calls),fixed=TRUE),logical(1))))
prov<-jsonlite::read_json(r$provenance,simplifyVector=TRUE)
stopifnot(setequal(unlist(prov$unsupported_accession_ids),ids[-1]),setequal(unlist(prov$requested_accession_ids),ids))
offline<-function(...)stop("Unexpected HTTP")
cached<-pm_reference_resolve("uniprot",cache_dir=file.path(td,"cache"),accessions=ids,fetch=offline)
stopifnot(cached$mode=="cache",identical(cached$snapshot,r$snapshot))
f<-data.frame(record_id=seq_along(ids),accession=ids,sequence="MACDEFGHIK")
m<-data.frame(record_id=seq_along(ids),accession=ids,peptide_id=paste0("p",seq_along(ids)),start=2L,end=5L)
a<-pm_annotate_uniprot(list(fasta=f,matches=m),r$snapshot,r$provenance)
stopifnot(identical(a$matches,m),a$uniprot$records$status[1]=="matched",
  all(a$uniprot$records$status[-1]=="unsupported_identifier"))
empty<-pm_reference_resolve("uniprot",cache_dir=file.path(td,"all_custom"),accessions=ids[-1],fetch=offline)
ep<-jsonlite::read_json(empty$provenance,simplifyVector=TRUE)
stopifnot(ep$records==0L,ep$retrieval_status=="no_queryable_accessions",
  !any(startsWith(basename(list.dirs(file.path(td,"all_custom","uniprot"),recursive=FALSE)),".stage-")))
# Full annotated run keeps every mapped record, consumes the mixed cache offline,
# and writes explicit unsupported statuses rather than a false negative match.
writeLines(as.vector(rbind(paste0(">",ids),rep("MACDEFGHIK",length(ids)))),file.path(td,"input.fasta"))
write.csv(data.frame(Sequence="ACDE","Master Protein Accessions"=ids,check.names=FALSE),file.path(td,"input.csv"),row.names=FALSE)
cfg<-list(peptides="input.csv",fasta="input.fasta",output_dir="results",sequence_column="Sequence",excel=TRUE,cache_dir="cache")
saved_fetch<-pm_reference_fetch;pm_reference_fetch<-offline
pm_run(cfg,td,file.path(root,"inst","templates","report_template.html"))
pm_reference_fetch<-saved_fetch
records<-pm_read_table(file.path(td,"results","audit","uniprot_records.csv"))
matches<-pm_read_table(file.path(td,"results","audit","matches.csv"))
manifest<-jsonlite::read_json(file.path(td,"results","run_manifest.json"),simplifyVector=TRUE)
stopifnot(nrow(matches)==length(ids),setequal(matches$accession,ids),
  sum(records$status=="unsupported_identifier")==length(ids)-1L,
  manifest$reference_access$status=="partial_identifier_coverage")
# Supplied snapshots remain authoritative: custom IDs are usable if an explicit
# exact-ID/full-sequence record is provided. No automatic ID rewriting is added.
jsonlite::write_json(list(results=list(list(primaryAccession="CON__TEST",sequence=list(value="MACDEFGHIK"),features=list()))),file.path(td,"custom.json"),auto_unbox=TRUE)
jsonlite::write_json(list(schema_version="proteomaster-uniprot-rest-1",records=1L,
  derived_sha256=digest::digest(file=file.path(td,"custom.json"),algo="sha256")),file.path(td,"custom_prov.json"),auto_unbox=TRUE)
supplied<-pm_annotate_uniprot(list(fasta=f,matches=m),file.path(td,"custom.json"),file.path(td,"custom_prov.json"))
stopifnot(supplied$uniprot$records$status[supplied$uniprot$records$accession=="CON__TEST"]=="matched")
html<-pm_report_reference_html(list(enabled=TRUE,providers=list(),unsupported_uniprot_ids=list("custom_<tag>")))
stopifnot(grepl("custom_&lt;tag&gt;",html,fixed=TRUE),!grepl("custom_<tag>",html,fixed=TRUE))
cat("Mixed/custom identifiers: scoped queries, explicit status, unchanged mappings and offline cache passed.\n")
