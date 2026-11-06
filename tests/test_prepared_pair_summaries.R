args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1]),winslash="/")))
source(file.path(root,"R/report.R"))
source(file.path(root,"R/excel_locale.R"));source(file.path(root,"R/excel.R"))
source(file.path(root,"R/export_bundle.R"))
d <- jsonlite::fromJSON(file.path(root,"tests/fixtures/export_bundle.json"))
tables <- pm_export_tables(d)
s <- attr(tables,"pair_summaries")
stopifnot(nrow(s)>0L,all(s$arithmetic_valid),all(s$selection_rank[s$Representative]==1L))
# A compact report can regenerate every workbook sheet without ratio rows.
compact <- d
compact$pair_view$prepared_summaries <- jsonlite::fromJSON(jsonlite::toJSON(s,dataframe="rows",auto_unbox=TRUE,na="null",null="null",digits=I(17)))
compact$pair_view$ratios <- NULL
compact$pair_view$contrasts <- NULL
refreshed <- pm_export_tables(compact)
for (name in names(tables)) stopifnot(isTRUE(all.equal(tables[[name]],refreshed[[name]],check.attributes=FALSE,tolerance=0)))
stopifnot(isTRUE(all.equal(attr(tables,"overview"),attr(refreshed,"overview"),check.attributes=FALSE,tolerance=0)))
# Counts are sufficient to preserve one-condition/zero/missing states when users
# change display coverage without recalculating or changing the representative.
stopifnot(any(grepl("only in",s$status)),all(s$control_common<=s$control_neo),
  all(s$control_common<=s$control_spanning),all(s$treatment_common<=s$treatment_neo),
  all(s$treatment_common<=s$treatment_spanning))
paircols <- c("site_id","role","neo_peptide_id","spanning_peptide_id")
links <- unique(d$pair_view$ratios[paircols])
links$assignment_eligible <- TRUE
links$assignment_eligible[1] <- FALSE
rejected <- d; rejected$pair_view$links <- links
rt <- pm_export_tables(rejected); rs <- attr(rt,"pair_summaries")
stopifnot(nrow(rs)==nrow(s),sum(!rs$assignment_eligible)==nrow(d$pair_view$comparisons),
  !any(rs$Representative[!rs$assignment_eligible]),all(is.na(rs$selection_rank[!rs$assignment_eligible])),
  nrow(rt$`All pairs`)==nrow(tables$`All pairs`)-nrow(d$pair_view$comparisons))
# Mapping duplicates must not conceal conflicting measurements.
bad <- d
duplicate <- bad$pair_view$ratios[1,,drop=FALSE]
duplicate$site_id <- "alternative:mapping"
duplicate$neo_value <- duplicate$neo_value*2
bad$pair_view$ratios <- rbind(bad$pair_view$ratios,duplicate)
stopifnot(inherits(try(pm_report_pair_summaries(bad),silent=TRUE),"try-error"))
missing_row <- d; missing_row$pair_view$ratios <- missing_row$pair_view$ratios[-1,,drop=FALSE]
stopifnot(inherits(try(pm_report_pair_summaries(missing_row),silent=TRUE),"try-error"))
bad_ratio <- d; bad_ratio$pair_view$ratios$neo_value[1] <- 0
stopifnot(inherits(try(pm_report_pair_summaries(bad_ratio),silent=TRUE),"try-error"))
bad <- compact; bad$pair_view$prepared_summaries$control_neoMissing <- NULL
stopifnot(inherits(try(pm_report_pair_summaries(bad),silent=TRUE),"try-error"))
no_contrast <- d; no_contrast$pair_view$comparisons <- data.frame()
nt <- pm_export_tables(no_contrast); ns <- attr(nt,"pair_summaries")
stopifnot(nrow(ns)>0L,all(is.na(ns$comparison)),all(ns$status=="No comparison"),
  all(ns$control_common==0L),all(ns$treatment_common==0L),all(ns$selection_k==0L),
  nrow(nt$Results)>0L,nrow(nt$`All pairs`)==nrow(ns),
  all(nt$Results$Status=="No comparison"),all(is.na(nt$Results[["Ratio log2FC"]])))
no_contrast$pair_view$prepared_summaries <- ns; no_contrast$pair_view$ratios <- NULL
nrt <- pm_export_tables(no_contrast)
for(name in names(nt))stopifnot(isTRUE(all.equal(nt[[name]],nrt[[name]],check.attributes=FALSE,tolerance=0)))
cat("Prepared summaries: all-sheet exact round trip, excluded links, missingness counts and malformed-input checks passed.\n")
