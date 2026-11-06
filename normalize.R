#!/usr/bin/env Rscript
# Explicit software-neutral tables -> sample-level ratios, without a contrast.
args <- commandArgs(trailingOnly=TRUE)
if (.Platform$OS.type=="windows") suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8"))
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),winslash="/",mustWork=TRUE))
for(module in c("mapping","report","excel_locale","input_adapter","maxquant_input","sample_ratios","normalization_io")) source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
if(length(args)!=2L || args[1]!="--config") stop("Usage: Rscript normalize.R --config normalization.json")
for(package in c("jsonlite","digest")) if(!requireNamespace(package,quietly=TRUE)) stop("Install ",package," first.")
config_path <- normalizePath(args[2],winslash="/",mustWork=TRUE)
config <- jsonlite::fromJSON(config_path,simplifyVector=TRUE)
required <- c("peptides","proteins","samples","output_dir")
if(!all(required %in% names(config))) stop("Config requires peptides, proteins, samples and output_dir.")
adapter <- if(is.null(config$adapter)) "generic" else config$adapter
if(!is.character(adapter)||length(adapter)!=1L||is.na(adapter)||!adapter %in% c("generic","maxquant_peptides")) stop("adapter must be generic or maxquant_peptides.")
if(adapter=="generic" && is.null(config$spec)) stop("Generic input requires spec.")
if(adapter=="maxquant_peptides" && !is.null(config$spec)) stop("Use maxquant_options instead of spec with maxquant_peptides.")
resolve <- function(path) {
  if(!is.character(path) || length(path)!=1L || is.na(path) || !nzchar(path)) stop("File paths must be single nonempty strings.")
  if(grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)",path)) path else file.path(dirname(config_path),path)
}
for(name in c("peptides","proteins","samples","output_dir")) config[[name]] <- resolve(config[[name]])
excel <- if(is.null(config$excel)) TRUE else config$excel
if(!is.logical(excel) || length(excel)!=1L || is.na(excel)) stop("excel must be true or false.")
if(excel && !requireNamespace("openxlsx",quietly=TRUE)) stop('Install openxlsx or set "excel": false.')
if(dir.exists(config$output_dir) && length(list.files(config$output_dir,all.files=TRUE,no..=TRUE))) stop("Choose an empty output directory.")
message("Reading explicitly selected columns...")
peptides <- pm_read_table(config$peptides,config$peptide_sheet)
proteins <- pm_read_table(config$proteins,config$protein_sheet)
sample_sheet <- if(is.null(config$sample_sheet) && grepl("\\.xlsx$",config$samples,ignore.case=TRUE)) "Samples" else config$sample_sheet
samples <- pm_read_table(config$samples,sample_sheet)
if(adapter=="generic") {
  if(!is.null(config$maxquant_options)) stop("maxquant_options requires adapter=maxquant_peptides.")
  data <- pm_import_quant_tables(peptides,proteins,config$spec,samples)
} else {
  options <- config$maxquant_options
  if(is.null(options)) options <- list()
  if(!is.list(options)||any(!names(options) %in% c("protein_quantity","exclude_only_identified_by_site"))) stop("Unknown maxquant_options; use protein_quantity and exclude_only_identified_by_site.")
  data <- do.call(pm_import_maxquant_tables,c(list(peptides=peptides,proteins=proteins,samples=samples),options))
}
message("Calculating per-sample peptide/protein ratios...")
pm_validate_quantity_columns(peptides, data$samples$peptide_column)
pm_validate_quantity_columns(proteins, data$samples$protein_column, "Protein table")
ratios <- pm_sample_ratios(data)
inputs <- list(peptides=config$peptides,proteins=config$proteins,samples=config$samples,configuration=config_path)
pm_export_sample_ratios(data,ratios,config$output_dir,inputs,excel=excel,notes=config$notes)
message("Done. Results: ",config$output_dir)
