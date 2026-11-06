script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))source(f,encoding="UTF-8")
fails <- function(expr,pattern){e<-tryCatch({force(expr);NULL},error=identity);stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e),fixed=TRUE))}
td <- tempfile("pipeline_columns_");pm_wizard_demo(td)
loaded <- jsonlite::read_json(file.path(td,"mapping.json"),simplifyVector=TRUE)
loaded$excel <- TRUE
input <- pm_read_table(file.path(td,"peptides.csv"))
# Include a literal identifier and untouched source metadata in the round trip.
input$Literal <- c("001","1E8","a;b","\u03b1","", "x,y")
input[["Modification detail"]][2] <- "Oxidation [M1]"
for(i in 1:6) {
  # Two neo forms vary independently of their measured spanning forms.
  input[[paste0("Signal ",i)]][2] <- as.character(c(40,55,45,95,120,90)[i])
  input[[paste0("Signal ",i)]][5] <- as.character(c(25,45,35,90,95,85)[i])
}
pm_write_csv(input,file.path(td,"peptides.csv"))
input <- pm_read_table(file.path(td,"peptides.csv"))
template <- file.path(root,"inst","templates","report_template.html")
export_real <- pm_export;export_observed <- NULL
pm_export <- function(result,...) {export_observed <<- result;export_real(result,...)}
manifest <- pm_run(loaded,config_base=td,template=template)
stopifnot(file.exists(file.path(td,"results","atelente_audit.zip")),
  file.exists(file.path(td,"results","atelente_results.xlsx")),
  identical(openxlsx::getSheetNames(file.path(td,"results","atelente_results.xlsx")),c("Results","Peptides","Settings")))
stopifnot(file.exists(file.path(td,"results","report.html")),manifest$peptide_rows==6L,
  identical(export_observed$source_input,input),
  identical(export_observed$column_mapping_audit$fields$sequence,"Peptide sequence"),
  identical(export_observed$settings$assignment_mode,"representative"),
  identical(names(export_observed$source_peptides)[1],"Sequence"))
exported <- pm_read_table(file.path(td,"results","audit","source_input.csv"))
expected_cells <- input
# pm_read_table records the file path being read; the export has its own path.
attr(exported,"source_file") <- NULL;attr(expected_cells,"source_file") <- NULL
stopifnot(identical(exported,expected_cells))
audit <- jsonlite::read_json(file.path(td,"results","audit","column_mapping_audit.json"),simplifyVector=TRUE)
stopifnot(identical(audit$fields$sequence,"Peptide sequence"),audit$row_count==nrow(input),
  identical(audit$quantities$source_column,paste0("Signal ",1:6)))
mapped_samples <- pm_read_table(file.path(td,"results","audit","neo_spanning_samples.csv"))
stopifnot(identical(mapped_samples$Sample,paste0("Sample_",1:6)),
  identical(mapped_samples[["Peptide column"]],sprintf("Abundance: Selected_%03d",1:6)))
no_compare <- loaded;no_compare$output_dir <- "no_comparison"
no_compare$sample_design$comparisons <- NULL
no_compare$sample_design$samples$Condition <- ""
pm_run(no_compare,config_base=td,template=template)
stopifnot(file.exists(file.path(td,"no_comparison","report.html")),
  !file.exists(file.path(td,"no_comparison","audit","neo_spanning_contrasts.csv")),
  is.null(export_observed$pair_statistics))
no_cmp_samples <- pm_read_table(file.path(td,"no_comparison","audit","neo_spanning_samples.csv"))
stopifnot(all(is.na(no_cmp_samples$Condition)|!nzchar(no_cmp_samples$Condition)),identical(no_cmp_samples$Sample,paste0("Sample_",1:6)))
bad_config <- loaded;bad_config$neo_spanning_samples <- "samples.xlsx"
fails(pm_run(bad_config,config_base=td,template=template),"not both")

# Optional inference must not disable the same representative-assignment gates.
# Shared master sets stay eligible, while a disjoint neo/spanning set is audited.
shared_input<-input
shared_input[["Representative ID"]]<-"TEST1;TEST2"
shared_input[["Reported IDs"]]<-"TEST1;TEST2"
shared_input[["Representative ID"]][5]<-"TEST1"
shared_input[["Representative ID"]][4]<-"TEST2"
pm_write_csv(shared_input,file.path(td,"shared.csv"))
writeLines(c(">sp|TEST1|SYNTHETIC","MACDEKFGHIRK",">sp|TEST2|SYNTHETIC","MACDEKFGHIRK"),file.path(td,"shared.fasta"))
descriptive<-loaded;descriptive$peptides<-"shared.csv";descriptive$fasta<-"shared.fasta"
descriptive$output_dir<-"descriptive_assignment"
descriptive$sample_design$samples[["Experimental unit"]]<-""
pm_run(descriptive,config_base=td,template=template)
scope<-export_observed$pair_statistics
stopifnot(identical(scope$settings$test,FALSE),nrow(scope$measurement_results)==0L,
  all(c("TEST1","TEST2") %in% scope$links$accession[scope$links$scope_eligible &
    scope$links$neo_peptide_id=="pep_0000002" & scope$links$spanning_peptide_id=="pep_0000003"]),
  any(scope$links$exclusion_reason=="No shared master accession" &
    scope$links$neo_peptide_id=="pep_0000005" & scope$links$spanning_peptide_id=="pep_0000004"),
  file.exists(file.path(td,"descriptive_assignment","audit","pair_links.csv")))
descriptive_links<-pm_read_table(file.path(td,"descriptive_assignment","audit","pair_links.csv"))
stopifnot(any(descriptive_links$exclusion_reason=="No shared master accession",na.rm=TRUE))
descriptive_html<-paste(readLines(file.path(td,"descriptive_assignment","report.html"),warn=FALSE,encoding="UTF-8"),collapse="\n")
marker<-'<script type="application/json" id="data">'
marker_start<-regexpr(marker,descriptive_html,fixed=TRUE)[1]
stopifnot(marker_start>0L)
payload_tail<-substring(descriptive_html,marker_start+nchar(marker))
payload_end<-regexpr("</script>",payload_tail,fixed=TRUE)[1]
stopifnot(payload_end>0L)
descriptive_payload<-pm_report_decode_json(substr(payload_tail,1,payload_end-1L))
view<-descriptive_payload$pair_view
stopifnot(is.data.frame(view$links),nrow(view$links)==nrow(scope$links),
  "assignment_eligible" %in% names(view$links),
  any(view$links$assignment_eligible & view$links$neo_peptide_id=="pep_0000002" &
    view$links$spanning_peptide_id=="pep_0000003"),
  all(!view$links$assignment_eligible[view$links$exclusion_reason=="No shared master accession"]),
  length(view$statistics)==0L,identical(view$inference_settings$test,FALSE))
unique_descriptive<-loaded;unique_descriptive$assignment_mode<-"unique"
unique_descriptive$output_dir<-"unique_descriptive"
pm_run(unique_descriptive,config_base=td,template=template)
stopifnot(is.null(export_observed$pair_statistics))

# The stub below verifies dispatch and serialization only, never model validity.
actual_workflow <- if(exists("pm_pair_workflow",mode="function"))pm_pair_workflow else NULL
dispatched <- NULL
pm_pair_workflow <- function(result,design,experimental_units,assignment_mode,representative_column,abundance_type,test) {
  dispatched <<- list(result=result,design=design,experimental_units=experimental_units,assignment_mode=assignment_mode,
    representative_column=representative_column,abundance_type=abundance_type,test=test)
  list(results=data.frame(measurement_pair_id=character(),family=character(),comparison=character(),
    status=character(),log2_effect=numeric(),p_value=numeric(),adjusted_p_value=numeric()),
    excluded=data.frame(reason=character()),settings=list(dispatch_test_only=TRUE),warnings=character())
}
with_test <- loaded;with_test$output_dir <- "test_dispatch";with_test$pair_statistics <- TRUE
pm_run(with_test,config_base=td,template=template)
stopifnot(identical(dispatched$assignment_mode,"representative"),isTRUE(dispatched$test),
  identical(dispatched$representative_column,"Master Protein Accessions"),identical(dispatched$abundance_type,"abundance"),
  identical(dispatched$experimental_units[["Experimental unit"]],paste0("Independent_",1:6)),
  identical(dispatched$design$comparisons$Treatment,"Treatment"),
  file.exists(file.path(td,"test_dispatch","audit","pair_results.csv")))
settings <- jsonlite::read_json(file.path(td,"test_dispatch","audit","pair_statistics_settings.json"))
stopifnot(isTRUE(settings$settings$dispatch_test_only))
if(is.null(actual_workflow)) {
  rm(pm_pair_workflow)
  cat("Actual statistical workflow unavailable: dispatch/export verified with explicit stub only.\n")
} else {
  pm_pair_workflow <- actual_workflow
  with_test$output_dir <- "actual_tests"
  pm_run(with_test,config_base=td,template=template)
  stopifnot(file.exists(file.path(td,"actual_tests","audit","pair_statistics_settings.json")))
  actual_results <- pm_read_table(file.path(td,"actual_tests","audit","pair_measurement_results.csv"))
  stopifnot(sum(actual_results$status=="tested")>=2L)
  cat("Actual statistical workflow completed through pm_run.\n")
}
pm_export <- export_real
cat("Column-mapped pipeline, source preservation, optional comparison and statistics dispatch tests passed.\n")
