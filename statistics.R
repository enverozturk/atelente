#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(.Platform$OS.type=="windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8")))
root <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),winslash="/",mustWork=TRUE))
for(module in c("mapping","report","excel_locale","statistics","statistics_io")) source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
if(length(args)!=2L||args[1]!="--config") stop("Usage: Rscript statistics.R --config statistics.json")
for(package in c("jsonlite","digest","limma")) if(!requireNamespace(package,quietly=TRUE)) stop("Install ",package," first.")
config_path <- normalizePath(args[2],winslash="/",mustWork=TRUE)
config <- jsonlite::fromJSON(config_path,simplifyVector=TRUE)
if(!all(c("ratios_dir","samples","comparisons","output_dir","design") %in% names(config))) stop("Require ratios_dir, samples, comparisons, output_dir and explicit design.")
resolve <- function(path) {
  if(!is.character(path)||length(path)!=1L||is.na(path)||!nzchar(path)) stop("Paths must be nonempty strings.")
  if(grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)",path)) path else file.path(dirname(config_path),path)
}
for(name in c("ratios_dir","samples","comparisons","output_dir")) config[[name]] <- resolve(config[[name]])
excel <- if(is.null(config$excel)) TRUE else config$excel
if(!is.logical(excel)||length(excel)!=1L||is.na(excel)) stop("excel must be true or false.")
if(excel&&!requireNamespace("openxlsx",quietly=TRUE)) stop("Install openxlsx or disable Excel export.")
if(dir.exists(config$output_dir)&&length(list.files(config$output_dir,all.files=TRUE,no..=TRUE))) stop("Choose an empty output directory.")
input <- pm_load_ratio_statistics(config$ratios_dir)
sheet <- function(path,configured,default) if(is.null(configured)&&grepl("\\.xlsx$",path,ignore.case=TRUE)) default else configured
samples <- pm_read_table(config$samples,sheet(config$samples,config$sample_sheet,"Samples"))
comparisons <- pm_read_table(config$comparisons,sheet(config$comparisons,config$comparison_sheet,"Comparisons"))
message("Fitting explicit independent-unit comparisons...")
result <- pm_test_ratios(input$ratios,samples,comparisons,design=config$design)
pm_export_statistics(result,config$output_dir,c(input$inputs,list(samples=config$samples,comparisons=config$comparisons,configuration=config_path)),excel=excel,notes=config$notes)
message("Done. Results: ",config$output_dir)
