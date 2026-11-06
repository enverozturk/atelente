args<-commandArgs();root<-dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1]),winslash="/")))
source(file.path(root,"R/report.R"))
source(file.path(root,"R/excel_locale.R"));source(file.path(root,"R/excel.R"));source(file.path(root,"R/export_bundle.R"))
fixture<-file.path(root,"tests/fixtures/export_bundle.json")
d<-jsonlite::fromJSON(fixture);t<-pm_export_tables(d);expected<-jsonlite::fromJSON(file.path(root,"tests/fixtures/export_bundle_expected.json"))
stopifnot(identical(names(t),c("Results","Peptides","Settings","All pairs","Other candidates","Annotations")))
actual<-t$Results;i<-match(expected$neo,actual$`Neo form ID`)
stopifnot(nrow(actual)==nrow(expected),!anyNA(i),identical(actual$`Spanning form ID`[i],expected$span),isTRUE(all.equal(actual$`Ratio log2FC`[i],expected$fc,tolerance=1e-10)))
stopifnot(!anyDuplicated(t$Peptides$`Form ID`),all(unique(d$cleavage$evidence$peptide_id) %in% t$Peptides$`Form ID`))
d_custom<-d
d_custom$pair_view$selection<-list(k_by_comparison=as.list(stats::setNames(rep(999L,nrow(d$pair_view$comparisons)),d$pair_view$comparisons$Name)))
custom<-pm_export_tables(d_custom)
stopifnot(nrow(custom$Results)>0L,!"Passes default coverage" %in% names(custom$Results))
out<-tempfile(fileext=".xlsx");pm_write_excel_tables(t,out,c("Results","Peptides","Settings"));stopifnot(identical(openxlsx::getSheetNames(out),c("Results","Peptides","Settings")))
cat("Export bundle: JS representative parity, complete unique forms and default sheets passed.\n")

# Reader schema is concise, typed, and traceable; offline exports keep an unfiltered catalog.
overview<-attr(t,"overview")
stopifnot(is.data.frame(overview),nrow(overview)>=nrow(t$Results),
  all(c("Cleavage window","Neo modifications","Spanning modifications","Pattern","Available pairs") %in% names(t$Results)),
  !any(c("Representative","Passes default coverage","Site IDs","Export key","Control","Treatment") %in% names(t$Results)),
  identical(tail(names(t$Results),2),c("Neo form ID","Spanning form ID")),
  "Control measured (of 5)" %in% names(t$Results),is.numeric(t$Results[["Control measured (of 5)"]]),
  all(grepl("^Neo-[NC]$",overview$Role)),all(overview$`Available pairs`>=1L),
  all(vapply(overview$`Export key`,function(k)length(jsonlite::fromJSON(k))==4L,logical(1))),
  all(t$Peptides$Roles %in% c("Neo-N","Neo-C","Spanning")),
  !any(grepl("^[{[]",t$Settings$Value)),
  all(c("Frozen selection k","Default display coverage k") %in% t$Settings$Setting))
# Numeric input zero, tiny values and missingness remain distinct in sample columns.
q<-d; q$pair_view$measurements$value[1:3]<-c(0,1.25e-12,NA_real_)
q$pair_view$measurements$value_status[1:3]<-c("observed","observed","missing")
qt<-pm_export_tables(q);id<-q$pair_view$measurements$peptide_id[1];row<-match(id,qt$Peptides$`Form ID`)
quantity_label<-function(i)paste0("Quantity: ",q$pair_view$samples$Sample[match(q$pair_view$measurements$source_column[i],q$pair_view$samples[["Peptide column"]])])
stopifnot(qt$Peptides[[quantity_label(1)]][row]==0,
  qt$Peptides[[quantity_label(2)]][match(q$pair_view$measurements$peptide_id[2],qt$Peptides$`Form ID`)]==1.25e-12,
  is.na(qt$Peptides[[quantity_label(3)]][row]))
# Exported comparison count controls coverage headers; multiple contrasts stay rectangular.
multi<-d;cp<-d$pair_view$comparisons;cp$Name<-"Second comparison";multi$pair_view$comparisons<-rbind(multi$pair_view$comparisons,cp)
mt<-pm_export_tables(multi)
stopifnot(all(c("Control measured","Treatment measured") %in% names(mt$Results)),
 !any(grepl("measured \\(of",names(mt$Results))),length(unique(mt$Results$Comparison))==2L)
# Shared mappings count the same alternative form pair once, not once per accession.
shared<-d;first<-d$cleavage$sites$site_id[1];newsite<-"duplicate:mapping"
ss<-d$cleavage$sites[d$cleavage$sites$site_id==first,,drop=FALSE];ss$site_id<-newsite;ss$accession<-"ALTERNATIVE";shared$cleavage$sites<-rbind(shared$cleavage$sites,ss)
ee<-d$cleavage$evidence[d$cleavage$evidence$site_id==first,,drop=FALSE];ee$site_id<-newsite;ee$accession<-"ALTERNATIVE";shared$cleavage$evidence<-rbind(shared$cleavage$evidence,ee)
rr<-d$pair_view$ratios[d$pair_view$ratios$site_id==first,,drop=FALSE];rr$site_id<-newsite;rr$accession<-"ALTERNATIVE";shared$pair_view$ratios<-rbind(shared$pair_view$ratios,rr)
st<-pm_export_tables(shared);j<-match(t$Results$`Neo form ID`,st$Results$`Neo form ID`)
stopifnot(identical(t$Results$`Available pairs`,st$Results$`Available pairs`[j]),nrow(st$Results)==nrow(t$Results))
cat("Reader-facing schema, sample quantities, dynamic numeric coverage, and deduplicated alternatives passed.\n")

# Source headings survive the adapter; Settings is one row per sample.
named<-d
named$column_mapping_audit<-list(quantities=data.frame(canonical_column=d$pair_view$samples[["Peptide column"]],source_column=paste("Original",seq_len(nrow(d$pair_view$samples))),input_scale="log2"))
named$pair_view$samples$Sample<-paste("User sample",seq_len(nrow(d$pair_view$samples)))
named$merops<-list(matches=data.frame(site_id=d$cleavage$sites$site_id[1],protease="Fixture"))
nt<-pm_export_tables(named);design<-nt$Settings[nt$Settings$Section=="Samples",,drop=FALSE]
stopifnot(nrow(design)==nrow(named$pair_view$samples),
 identical(design$`Source column`,named$column_mapping_audit$quantities$source_column),
 all(design$`Input scale`=="log2"),all(design$`Export scale`=="linear"),
 identical(nt$Peptides[["Quantity: User sample 1"]],t$Peptides[[paste0("Quantity: ",d$pair_view$samples$Sample[1])]]),
 all(paste0("Quantity: ",named$pair_view$samples$Sample) %in% names(nt$Peptides)),
 identical(nt$Annotations$`Annotation ID`,nt$Annotations$site_id),
 nt$Annotations$`Annotation ID`[1] %in% unlist(strsplit(nt$Results$`Annotation IDs`,"; ",fixed=TRUE)))
cat("Original/sample headings, rectangular sample design and annotation join keys passed.\n")

# Wholly missing components are audit-only, not reader results or Other candidates.
# The existing fixture also covers neo absent in control but jointly measured
# in treatment, and both components measured only in treatment: retain both.
stopifnot(!"span_only_neo" %in% t$Results$`Neo form ID`,
  !"span_only_neo" %in% t$`Other candidates`$`Form ID`,
  "span_only_neo" %in% t$`All pairs`$`Neo form ID`,
  all(c("neo_only_neo","both_only_neo") %in% t$Results$`Neo form ID`))
reverse<-d
ii<-which(reverse$pair_view$ratios$neo_peptide_id=="span_only_neo")
for(field in c("value","status")){
  a<-paste0("neo_",field);b<-paste0("spanning_",field)
  old<-reverse$pair_view$ratios[[a]][ii]
  reverse$pair_view$ratios[[a]][ii]<-reverse$pair_view$ratios[[b]][ii]
  reverse$pair_view$ratios[[b]][ii]<-old
}
mm<-reverse$pair_view$measurements
mm$value[mm$peptide_id=="span_only_neo"]<-10
mm$value_status[mm$peptide_id=="span_only_neo"]<-"observed"
mm$value[mm$peptide_id=="span_only_span"]<-NA_real_
mm$value_status[mm$peptide_id=="span_only_span"]<-"missing"
reverse$pair_view$measurements<-mm
rt<-pm_export_tables(reverse)
stopifnot(!"span_only_neo" %in% rt$Results$`Neo form ID`,
  !"span_only_neo" %in% rt$`Other candidates`$`Form ID`,
  nrow(rt$`All pairs`)==nrow(t$`All pairs`),
  "Neo quantified; spanning not quantified" %in% rt$`All pairs`$Status)
cat("Missing-component rows remain audit-only; one-condition results preserved.\n")
