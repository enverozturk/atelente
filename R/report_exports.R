# Reader-facing output bundle; the audit remains a separate on-disk ZIP.
pm_prepare_export_bundle <- function(payload, output_dir, manifest=list(), excel=TRUE) {
  payload$pair_view <- pm_report_selection(payload$pair_view)
  payload$export_settings <- manifest[c("settings","metadata","inputs","version","report_version","reference_access","pair_statistics")]
  if (isTRUE(excel)) {
    tables <- pm_export_tables(payload)
    if (!is.null(payload$pair_view)) payload$pair_view$prepared_summaries <- attr(tables,"pair_summaries")
    pm_write_excel_tables(tables,file.path(output_dir,"atelente_results.xlsx"),c("Results","Peptides","Settings"))
    dir.create(file.path(output_dir,"audit"),showWarnings=FALSE)
    all_path <- file.path(output_dir,"audit","atelente_tables.xlsx")
    public_sheets <- names(tables)
    overview <- attr(tables,"overview")
    internal_sheets <- character()
    if(is.data.frame(overview)) {
      tables[["__Results"]] <- overview
      internal_sheets <- "__Results"
    }
    pm_write_excel_tables(tables,all_path,names(tables))
    encoded <- gsub("[\r\n]","",jsonlite::base64_enc(readBin(all_path,"raw",n=file.info(all_path)$size)))
    payload$export_workbook <- list(sheets=public_sheets,internal_sheets=I(internal_sheets),base64=encoded)
    payload$excel_file <- "atelente_results.xlsx"
  } else {
    payload$excel_file <- NULL; payload$export_workbook <- NULL
  }
  payload$audit_file <- "atelente_audit.zip"
  payload$sample_metadata_file <- NULL
  # Audit references are retained for traceability and backward-compatible refresh.
  prefix <- function(path) {
    if (!is.character(path)||length(path)!=1L||is.na(path)||!nzchar(path)) return(path)
    if(startsWith(path,"audit/")) path else paste0("audit/",path)
  }
  if(length(payload$analysis_files)) payload$analysis_files <- lapply(payload$analysis_files,prefix)
  if(length(payload$quantity_files)) payload$quantity_files <- lapply(payload$quantity_files,prefix)
  if(!is.null(payload$evidence_workbook)) payload$evidence_workbook <- prefix(payload$evidence_workbook)
  payload
}
