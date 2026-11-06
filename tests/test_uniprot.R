args<-commandArgs();root<-dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1]),winslash="/")))
source(file.path(root,"R","uniprot.R"))
td<-tempfile("uniprot_");dir.create(td);snapshot<-file.path(td,"snapshot.json.gz");prov<-file.path(td,"provenance.json")
seq1<-"MACDEFGHIKLNPQRSTVWY"
feature<-function(type,a,b,description="",am="EXACT",bm="EXACT",evidence=NULL,context=NULL) {
  loc<-list(start=list(value=a,modifier=am),end=list(value=b,modifier=bm))
  if(!is.null(context)) loc$sequence<-context
  list(type=type,location=loc,description=description,evidences=evidence)
}
experimental<-list(list(evidenceCode="ECO:0000269",source="PubMed",id="123"))
prediction<-list(list(evidenceCode="ECO:0000255",source="PROSITE-ProRule",id="test"))
features<-list(feature("Chain",2,20,"Mature chain"),feature("Signal",1,4,evidence=prediction),
  feature("Transit peptide",1,3),feature("Propeptide",6,8,"Internal propeptide"),
  feature("Peptide",9,12,"Released peptide"),feature("Initiator methionine",1,1,"Removed",evidence=experimental),
  feature("Initiator methionine",1,1,"Not removed"),feature("Signal",1,7,bm="UNKNOWN"),
  feature("Chain",9,NA,bm="UNKNOWN"),feature("Peptide",14,16,context="OTHER-2"),
  feature("Propeptide",12,10),feature("Chain",1,999),
  feature("Peptide",14,18,am="OUTSIDE",bm="UNSURE"),feature("Initiator methionine",1,1,""))
entry<-function(acc,sequence=seq1,fs=features) list(primaryAccession=acc,entryAudit=list(entryVersion=3,sequenceVersion=2,
  lastSequenceUpdateDate="2020-01-01",lastAnnotationUpdateDate="2025-01-01"),sequence=list(value=sequence),features=fs)
main<-entry("T1");main$secondaryAccessions<-list("OLD")
entries<-list(main,entry("DIFF",paste0(seq1,"A")),entry("DUP"),entry("DUP"),entry("ISO-2"),entry("EMPTY",fs=list()),entry("NOSEQ",sequence=""))
save_snapshot<-function() {
  text<-pm_up_json(list(results=entries));con<-gzfile(snapshot,"wb");writeBin(charToRaw(enc2utf8(text)),con);close(con)
  jsonlite::write_json(list(schema_version="proteomaster-uniprot-rest-1",derived_sha256=digest::digest(file=snapshot,algo="sha256"),records=length(entries)),prov,auto_unbox=TRUE)
}
save_snapshot()
fa<-data.frame(record_id=1:9,accession=c("T1","T1","DIFF","T1-2","NONE","DUP","ISO-2","OLD","EMPTY"),sequence=seq1,stringsAsFactors=FALSE)
m<-data.frame(peptide_id=paste0("p",1:12),record_id=c(1,1,1,1,2:9),start=c(5,1,2,9,rep(5,8)),end=c(8,20,12,10,rep(8,8)),stringsAsFactors=FALSE)
sites<-data.frame(site_id=paste0("s",1:12),record_id=c(rep(1,8),2,3,4,7),bond_left=c(1,3,4,5,7,8,12,15,4,4,4,4),stringsAsFactors=FALSE)
input<-list(fasta=fa,matches=m,cleavage=list(sites=sites,profile=list(test=TRUE)),measurements=data.frame(peptide_id=c("p1","p2"),value=c(0,NA_real_)))
out<-pm_annotate_uniprot(input,snapshot,prov);u<-out$uniprot
stopifnot(identical(out[names(input)],input),
  identical(u$records$status,c("matched","matched","sequence_mismatch","unsupported_isoform","not_found","ambiguous_accession","matched","not_found","matched")),
  all(u$records$entry_version[u$records$status=="matched"]=="3"),
  !any(u$site_annotations$record_id %in% c(3,4,5,6,8,9)),
  !"s5" %in% u$site_annotations$site_id, !"s8" %in% u$site_annotations$site_id,
  "s9" %in% u$site_annotations$site_id, "s12" %in% u$site_annotations$site_id,
  sum(u$site_annotations$site_id=="s1")==2L,
  sum(u$site_annotations$site_id=="s6")==3L,
  all(u$site_annotations$internal_bond),
  all(c("sequence_mismatch","ambiguous_accession","different_sequence_context","invalid_range","out_of_range","partial_exact_boundaries","uncertain_boundaries") %in% u$features$coordinate_status),
  "ECO:0000255" %in% u$site_annotations$evidence_codes,
  "ECO:0000269" %in% u$site_annotations$evidence_codes,
  !any(u$site_annotations$description=="Not removed"),
  !any(u$site_annotations$feature_type=="Initiator methionine" & u$site_annotations$description==""))
t<-u$terminal_annotations
stopifnot(any(t$peptide_id=="p1" & t$role=="N" & t$boundary=="after_signal"),
  any(t$peptide_id=="p1" & t$role=="C" & t$boundary=="after_propeptide"),
  any(t$peptide_id=="p2" & t$role=="C" & t$boundary=="chain_end" & !t$internal_bond),
  sum(t$peptide_id=="p3" & t$role=="N")==2L,
  any(t$peptide_id=="p3" & t$role=="C" & t$boundary=="peptide_end"))
# Full-length Chain contributes reference-end annotations but no artificial internal site.
natural<-entry("NAT",fs=list(feature("Chain",1,20,"Full length")));entries<-list(natural);save_snapshot()
nat<-list(fasta=data.frame(record_id=1,accession="NAT",sequence=seq1),matches=data.frame(record_id=1,peptide_id="p",start=1,end=20))
n<-pm_annotate_uniprot(nat,snapshot,prov)
stopifnot(nrow(n$uniprot$site_annotations)==0,nrow(n$uniprot$terminal_annotations)==2,
  all(!n$uniprot$terminal_annotations$internal_bond),identical(n[names(nat)],nat))
# No mapped observations produces empty typed tables, not artificial negatives.
nat$matches<-nat$matches[FALSE,];empty<-pm_annotate_uniprot(nat,snapshot,prov)
stopifnot(nrow(empty$uniprot$records)==0,nrow(empty$uniprot$features)==0,nrow(empty$uniprot$terminal_annotations)==0)
# Observed alternate Met removal is retained with its qualifier; isoform-specific
# locations cannot leak into an otherwise sequence-matched canonical record.
entries<-list(entry("ALT",fs=list(feature("Initiator methionine",1,1,"Removed; alternate"),
  feature("Initiator methionine",1,1,"Removed",context="ALT-5"),
  feature("Transit peptide",1,NULL,bm="UNKNOWN"))))
save_snapshot();alt<-list(fasta=data.frame(record_id=1,accession="ALT",sequence=seq1),matches=data.frame(record_id=1,peptide_id="alt",start=2,end=10))
ao<-pm_annotate_uniprot(alt,snapshot,prov)
stopifnot(nrow(ao$uniprot$terminal_annotations)==1,
  ao$uniprot$terminal_annotations$description=="Removed; alternate",
  ao$uniprot$features$coordinate_status[2]=="different_sequence_context",
  ao$uniprot$features$coordinate_status[3]=="partial_exact_boundaries")
# Tamper and structural mismatch fail before annotation.
cat("tampered",file=snapshot,append=TRUE)
stopifnot(inherits(try(pm_annotate_uniprot(nat,snapshot,prov),silent=TRUE),"try-error"))
save_snapshot();p<-jsonlite::fromJSON(prov);p$records<-99;jsonlite::write_json(p,prov,auto_unbox=TRUE)
stopifnot(inherits(try(pm_annotate_uniprot(nat,snapshot,prov),silent=TRUE),"try-error"))
cat("UniProt tests passed.\n")
