args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1]),winslash="/")))
source(file.path(root,"R","merops.R"))
td <- tempfile("merops_");dir.create(td)
snapshot <- file.path(td,"snapshot.csv.gz");prov <- file.path(td,"provenance.json")
reference <- data.frame(merops_record_id=paste0("CLE",1:6),protease_id=rep("C00.001",6),
  protease_name=rep("Synthetic test enzyme",6),substrate_accession=c("P1","P1","P2","P3","P1-2","P1"),
  p1_position=c("4","4","10","9","4","4"),substrate_species=rep("Test species",6),
  p4_p4prime=c("ACDEFGHI","ACDEFGHV","ACDEFGHI","ACDEFGHI","ACDEFGHI",""),
  canonical_octamer=c(rep("true",5),"false"),
  cleavage_type=c("physiological","non-physiological","synthetic","theoretical","pathological","physiological"),
  method_code_raw=c("MS","NT","","CS","",""),additional_evidence_code_raw=rep("",6),
  reference_raw=c("Test reference\nwith embedded tab\tretained \u03b1 \u0094",rep("",5)),reference_ids=rep("test",6),stringsAsFactors=FALSE)
save_snapshot <- function() {
  # Byte-level fixture writing avoids write.csv's Windows native-codepage
  # substitution; production snapshots are generated as UTF8 by the preparer.
  quote_field <- function(x) paste0('"',gsub('"','""',enc2utf8(x),fixed=TRUE),'"')
  lines <- c(paste(quote_field(names(reference)),collapse=","),vapply(seq_len(nrow(reference)),
    function(i) paste(quote_field(unlist(reference[i,],use.names=FALSE)),collapse=","),character(1)))
  con <- gzfile(snapshot,"wb");writeBin(charToRaw(enc2utf8(paste0(paste(lines,collapse="\n"),"\n"))),con);close(con)
  jsonlite::write_json(list(schema_version="proteomaster-merops-observed-export-1",
    derived_sha256=digest::digest(file=snapshot,algo="sha256"),records=nrow(reference)),prov,auto_unbox=TRUE)
}
save_snapshot()
sites <- data.frame(site_id=c("a","b","c","d","e"),accession=c("P1","P1-3","P1-2","NONE","P1"),
  bond_left=c(4,4,4,12,4),bond_right=c(5,5,5,13,5),
  context_display=c("ACDE|FGHI","ACDE|FGHI","ACDE|FGHI","YYYY|YYYY","---E|FGHI"),
  context_complete=c(TRUE,TRUE,TRUE,TRUE,FALSE),stringsAsFactors=FALSE)
input <- list(cleavage=list(sites=sites),measurements=data.frame(value=c(NA,0,10)),peptides=data.frame(id=1:3))
out <- pm_annotate_merops(input,snapshot,prov);m <- out$merops$matches
stopifnot(identical(input$measurements,out$measurements),identical(input$peptides,out$peptides),
  identical(input$cleavage,out$cleavage),nrow(m)==17L,
  m$context_status[m$site_id=="a" & m$merops_record_id=="CLE1"]=="agreement",
  m$context_status[m$site_id=="a" & m$merops_record_id=="CLE2"]=="mismatch",
  m$context_status[m$site_id=="a" & m$merops_record_id=="CLE6"]=="unknown",
  all(m$relation_type[m$site_id=="b"]=="motif_compatibility"),
  sum(m$relation_type[m$site_id=="c"]=="reference_site_candidate")==1L,
  !"d" %in% m$site_id,all(m$context_status[m$site_id=="e"]=="unknown"),
  all(m$relation_type[m$site_id=="e"]=="reference_site_candidate"),
  all(c("synthetic","theoretical","pathological") %in% m$cleavage_type),
  m$reference_raw[m$site_id=="a" & m$merops_record_id=="CLE1"]==reference$reference_raw[1])
# Empty sites retain a typed, exportable table; no artificial no-hit rows.
empty <- input;empty$cleavage$sites <- sites[FALSE,]
stopifnot(nrow(pm_annotate_merops(empty,snapshot,prov)$merops$matches)==0L)
# Mutated bytes must fail before being parsed.
cat("tamper",file=snapshot,append=TRUE)
stopifnot(inherits(try(pm_annotate_merops(input,snapshot,prov),silent=TRUE),"try-error"))
save_snapshot()
reference$canonical_octamer[1] <- "false";save_snapshot()
stopifnot(inherits(try(pm_annotate_merops(input,snapshot,prov),silent=TRUE),"try-error"))
reference$canonical_octamer[1] <- "true";reference$merops_record_id[2] <- "CLE1";save_snapshot()
stopifnot(inherits(try(pm_annotate_merops(input,snapshot,prov),silent=TRUE),"try-error"))
cat("MEROPS tests passed.\n")
