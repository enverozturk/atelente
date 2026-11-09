# Synthetic presentation metadata only: no analysis run or workbook generation.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
e <- new.env(parent=globalenv())
for(f in c("report.R","excel_locale.R","excel.R","export_bundle.R","report_exports.R"))sys.source(file.path(root,"R",f),envir=e)
evidence <- data.frame(site_id="focal",peptide_id="span",record_id=1L,accession="SAME_ACCESSION",start=15L,end=25L,role="spanning")
terms <- evidence[rep(1,2),setdiff(names(evidence),"site_id"),drop=FALSE]
terms$role <- c("N","C");terms$processing_status <- "not_assessed";terms$classification <- "digestion_compatible"
terms$pre_digest_n_label <- NA;terms$digestion_compatible <- TRUE
payload <- function(t=terms,ev=evidence)list(cleavage=list(evidence=ev,termini=t))
expect <- function(d,flag,text=NULL){
  x <- e$pm_report_spanning_ends(d)
  stopifnot(nrow(x)==1L,identical(x$terminal,flag))
  if(!is.null(text))stopifnot(grepl(text,x$ends,fixed=TRUE))
  invisible(x)
}
expect(payload(),"No","N: digestion-compatible; C: digestion-compatible")
for(natural in c("reference_end","possible_initiator_met_removal")){
  t <- terms;t$processing_status[1] <- natural;t$pre_digest_n_label[1] <- TRUE;t$digestion_compatible[1] <- FALSE
  expect(payload(t),"No",if(natural=="reference_end")"N: protein end" else "N: possible Met removal")
}
t <- terms;t$pre_digest_n_label[1] <- TRUE;expect(payload(t),"Yes","N: pre-digestion N label")
t <- terms;t$digestion_compatible[2] <- FALSE;expect(payload(t),"Yes","C: not digestion-compatible")
expect(payload(terms[1,,drop=FALSE]),"Unknown","C: unknown")
t <- terms[1,,drop=FALSE];t$digestion_compatible <- FALSE;expect(payload(t),"Yes","C: unknown")
expect(payload(terms[FALSE,,drop=FALSE]),"Unknown","N: unknown")
for(field in c("processing_status","classification","pre_digest_n_label")){
  t <- terms;t$digestion_compatible[1] <- FALSE;t[[field]] <- NULL
  expect(payload(t),"Yes","N: not digestion-compatible")
}
t <- terms;t$digestion_compatible <- NULL;expect(payload(t),"Unknown","N: unknown")
t <- terms;t$processing_status <- NULL;t$classification[1] <- "reference_end";expect(payload(t),"No","N: protein end")
t <- terms;t$classification <- NULL;t$processing_status[1] <- "reference_end";expect(payload(t),"No","N: protein end")
# Ignore a different placement unless it appears in focal spanning evidence.
tt <- terms;tt$start <- 41L;tt$end <- 51L;tt$digestion_compatible <- FALSE
expect(payload(rbind(terms,tt)),"No")
ev <- evidence;ev$start <- 14L;ev$end <- 24L
tt <- terms;tt$start <- 14L;tt$end <- 24L;tt$digestion_compatible[1] <- FALSE
x <- expect(payload(rbind(terms,tt),rbind(evidence,ev)),"Yes")
stopifnot(grepl("SAME_ACCESSION [1 15-25]",x$ends,fixed=TRUE),grepl("SAME_ACCESSION [1 14-24]",x$ends,fixed=TRUE))
stopifnot(identical(x,e$pm_report_spanning_ends(payload(rbind(terms,tt),rbind(evidence,ev,ev)))))
# All old fields remain identical when descriptive metadata changes. Freeze the
# already-computed representative summaries to isolate presentation behavior.
d <- jsonlite::fromJSON(file.path(root,"tests/fixtures/export_bundle.json"))
fixed <- e$pm_report_pair_summaries(d)
native <- fixed
e$pm_report_pair_summaries <- function(payload)fixed
before <- e$pm_export_tables(d)
d$cleavage$termini$pre_digest_n_label <- rep(TRUE,nrow(d$cleavage$termini))
after <- e$pm_export_tables(d)
added <- c("Spanning ends","Spanning terminal form")
for(sheet in names(before)){
  old <- setdiff(names(before[[sheet]]),added)
  stopifnot(identical(before[[sheet]][old],after[[sheet]][old]))
}
stopifnot(identical(attr(before,"pair_summaries"),attr(after,"pair_summaries")))
stopifnot(all(added %in% names(after$Results)),all(added %in% names(after$`All pairs`)))
# Neo-C is a focal role, not a restriction on which spanning end is described.
fixed$role <- "C"
c_results <- e$pm_export_tables(d)
stopifnot(all(c_results$`All pairs`$Role=="Neo-C"),all(c_results$`All pairs`$`Spanning terminal form`=="Yes"))
# The canonical HTML lookup is attached even when workbook generation is disabled.
canonical <- e$pm_report_spanning_ends(d)
prepared <- e$pm_prepare_export_bundle(d,tempdir(),excel=FALSE)
stopifnot(identical(prepared$pair_view$spanning_end_reporting,canonical),is.null(prepared$export_workbook))
# The same observation pair can map to distinct records with the same accession.
g <- jsonlite::fromJSON(file.path(root,"tests/fixtures/export_bundle.json"))
site <- native$site_id[1];span <- native$spanning_peptide_id[1];neo <- native$neo_peptide_id[1]
ii <- g$cleavage$termini$record_id==1L & g$cleavage$termini$peptide_id==span
g$cleavage$termini$digestion_compatible[ii] <- TRUE;g$cleavage$termini$pre_digest_n_label <- NA
ss <- g$cleavage$sites[g$cleavage$sites$site_id==site,,drop=FALSE];ss$site_id <- "alternate_record";ss$record_id <- 99L
g$cleavage$sites <- rbind(g$cleavage$sites,ss)
ee <- g$cleavage$evidence[g$cleavage$evidence$site_id==site,,drop=FALSE];ee$site_id <- ss$site_id;ee$record_id <- 99L
g$cleavage$evidence <- rbind(g$cleavage$evidence,ee)
tt <- g$cleavage$termini[g$cleavage$termini$record_id==1L,,drop=FALSE];tt$record_id <- 99L
tt$digestion_compatible[tt$peptide_id==span & tt$role=="N"] <- FALSE
g$cleavage$termini <- rbind(g$cleavage$termini,tt)
pp <- native[native$site_id==site,,drop=FALSE];pp$site_id <- ss$site_id;pp$record_id <- 99L
fixed <- rbind(native,pp)
gt <- e$pm_export_tables(g);lookup <- e$pm_report_spanning_ends(g)
stopifnot(lookup$terminal[lookup$site_id==site & lookup$spanning_peptide_id==span]=="No",
  lookup$terminal[lookup$site_id==ss$site_id & lookup$spanning_peptide_id==span]=="Yes")
rr <- gt$Results[gt$Results$`Neo form ID`==neo,,drop=FALSE]
stopifnot(nrow(rr)==1L,rr$`Spanning terminal form`=="Yes",grepl("[1 ",rr$`Spanning ends`,fixed=TRUE),grepl("[99 ",rr$`Spanning ends`,fixed=TRUE))
cat("Spanning metadata: precedence, missing fields/termini, focal and repeated placements, Neo-C, canonical lookup and unchanged old fields passed.\n")
