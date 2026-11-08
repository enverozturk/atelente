#!/usr/bin/env Rscript
# Add verified quantities to an existing mapping report; never change its source directory.
if(.Platform$OS.type=="windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8")))
args <- commandArgs(trailingOnly=TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),winslash="/",mustWork=TRUE)))
for(module in c("mapping","report","statistics_io","report_quantities")) source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
if(length(args)!=2L||args[1]!="--config") stop("Usage: Rscript scripts/refresh_report.R --config report.json")
for(package in c("jsonlite","digest")) if(!requireNamespace(package,quietly=TRUE)) stop("Install ",package," first.")
config_path <- normalizePath(args[2],winslash="/",mustWork=TRUE)
config <- jsonlite::fromJSON(config_path,simplifyVector=TRUE)
if(!all(c("mapping_dir","ratios_dir","output_dir") %in% names(config))) stop("Require mapping_dir, ratios_dir and output_dir.")
resolve <- function(path) {
  if(!is.character(path)||length(path)!=1L||is.na(path)||!nzchar(path)) stop("Paths must be nonempty strings.")
  if(grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)",path)) path else file.path(dirname(config_path),path)
}
for(name in intersect(c("mapping_dir","ratios_dir","statistics_dir","output_dir"),names(config))) config[[name]] <- resolve(config[[name]])
if(file.exists(config$output_dir)&&!dir.exists(config$output_dir)) stop("output_dir must be a directory.")
if(dir.exists(config$output_dir)&&length(list.files(config$output_dir,all.files=TRUE,no..=TRUE))) stop("Choose an empty output directory.")
mapping <- normalizePath(config$mapping_dir,winslash="/",mustWork=TRUE)
html_path <- file.path(mapping,"report.html")
html <- paste(readLines(html_path,encoding="UTF-8",warn=FALSE),collapse="\n")
marker <- '<script type="application/json" id="data">'
starts <- gregexpr(marker,html,fixed=TRUE)[[1]]
if(length(starts)!=1L||starts[1]<0L) stop("Mapping report must have one embedded data block.")
tail <- substring(html,starts[1]+nchar(marker),nchar(html)); end <- regexpr('</script>',tail,fixed=TRUE)[1]
if(end<1L) stop("Mapping data block is incomplete.")
payload <- pm_report_decode_json(substr(tail,1,end-1L),max_uncompressed_bytes=450*1024^2)
if(!is.null(payload$quantities)) stop("Use the original mapping report, without an attached quantity layer.")
artifact_root <- function(directory,required_file) {
  # Both historical flat exports and the compact audit/ layout are supported.
  if(file.exists(file.path(directory,required_file))) return(directory)
  candidate <- file.path(directory,"audit")
  if(file.exists(file.path(candidate,required_file))) {
    parent <- paste0(normalizePath(directory,winslash="/",mustWork=TRUE),"/")
    actual <- normalizePath(candidate,winslash="/",mustWork=TRUE)
    if(!startsWith(paste0(actual,"/"),parent)) stop("Audit path escapes its result directory.")
    return(candidate)
  }
  stop("Missing report artifact: ",required_file)
}
source_path <- file.path(artifact_root(mapping,"source_peptides.csv"),"source_peptides.csv")
original <- pm_read_table(source_path)
if(!identical(names(original)[1:2],c("peptide_id","source_row"))) stop("Expected original mapping source export with peptide_id and source_row prefix.")
if(!identical(as.character(original$peptide_id),as.character(payload$peptides$peptide_id))||
   !identical(as.character(original$source_row),as.character(payload$peptides$source_row))) stop("Mapping report/source observation order differs.")
result <- list(peptides=payload$peptides,source_peptides=original[,-c(1L,2L),drop=FALSE])
if ("Leading razor protein" %in% names(result$source_peptides))
  payload$peptides$leading_razor_accessions <- result$source_peptides[["Leading razor protein"]]
message("Verifying observation identity and quantity exports...")
ratios_root <- artifact_root(config$ratios_dir,"normalization_manifest.json")
statistics_root <- if(!is.null(config$statistics_dir)) artifact_root(config$statistics_dir,"statistics_manifest.json") else NULL
quantities <- pm_report_quantities(result,ratios_root,statistics_root)
# Copy adjacent generated artifacts unchanged; no recomputation or external fetch.
copy_artifacts <- function(from,to) {
  dir.create(to,recursive=TRUE,showWarnings=FALSE)
  paths <- list.files(from,pattern="\\.(csv|json|xlsx|zip|md)$",full.names=TRUE,ignore.case=TRUE,recursive=TRUE)
  paths <- paths[!dir.exists(paths)]
  base <- paste0(normalizePath(from,winslash="/",mustWork=TRUE),"/")
  actual <- normalizePath(paths,winslash="/",mustWork=TRUE)
  if(any(!startsWith(actual,base))) stop("Report artifact escapes its result directory.")
  relative <- substring(gsub("\\\\","/",paths),nchar(gsub("\\\\","/",from))+2L)
  for(i in seq_along(paths)) {
    path <- paths[i]
    dest <- file.path(to,relative[i])
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(path,dest,overwrite=FALSE)) stop("Could not copy ",basename(path))
    if(!identical(digest::digest(file=path,algo="sha256"),digest::digest(file=dest,algo="sha256"))) stop("Copy hash mismatch.")
  }
  stats::setNames(lapply(paths,function(p)digest::digest(file=p,algo="sha256")),relative)
}
mapping_hashes <- copy_artifacts(mapping,config$output_dir)
ratio_hashes <- copy_artifacts(ratios_root,file.path(config$output_dir,"quantities"))
statistics_hashes <- if(!is.null(statistics_root)) copy_artifacts(statistics_root,file.path(config$output_dir,"statistics")) else NULL
payload$quantities <- quantities
payload$quantity_files <- list(ratios="quantities/sample_ratios.csv",manifest="report_manifest.json")
if(file.exists(file.path(config$output_dir,"quantities","sample_ratios.xlsx"))) payload$quantity_files$workbook <- "quantities/sample_ratios.xlsx"
if(!is.null(config$statistics_dir)) {
  payload$quantity_files$statistics <- "statistics/results.csv"
  if(file.exists(file.path(config$output_dir,"statistics","peptide_statistics.xlsx"))) payload$quantity_files$statistics_workbook <- "statistics/peptide_statistics.xlsx"
  if(!is.null(payload$cleavage)) {
    paired <- pm_paired_bh(quantities$results,payload$cleavage)
    payload$quantities$results <- paired$results
    payload$quantities$paired_bh_counts <- paired$counts
    paired_file <- "paired_peptide_statistics.csv"
    paired_rows <- paired$results[paired$results$paired_evidence,c("peptide_id","observation_id","entity_id","comparison","sequence","status","log2_effect","p_value","paired_adjusted_p_value"),drop=FALSE]
    pm_write_csv(paired_rows,file.path(config$output_dir,paired_file))
    payload$quantity_files$paired_statistics <- paired_file
    quantities$provenance$paired_bh <- list(method="BH",family="tested observations supporting both-termini sites, once per comparison",counts=paired$counts,
      sha256=digest::digest(file=file.path(config$output_dir,paired_file),algo="sha256"))
  }
}
pm_write_report(payload,file.path(root,"inst","templates","report_template.html"),file.path(config$output_dir,"report.html"))
manifest <- list(version="0.8.0",stage="integrated_quantity_report",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
  mapping_report_sha256=digest::digest(file=html_path,algo="sha256"),mapping_source_sha256=digest::digest(file=source_path,algo="sha256"),
  provenance=quantities$provenance,unchanged_artifact_hashes=list(mapping=mapping_hashes,quantities=ratio_hashes,statistics=statistics_hashes),
  report_sha256=digest::digest(file=file.path(config$output_dir,"report.html"),algo="sha256"),
  scope="Observation-level presentation only; no statistical refit, mapping changes, site aggregation or protease attribution.")
if (!is.null(payload$reference_access)) manifest$reference_access <- payload$reference_access
jsonlite::write_json(manifest,file.path(config$output_dir,"report_manifest.json"),pretty=TRUE,auto_unbox=TRUE,na="null",null="null")
message("Done. Report: ",file.path(config$output_dir,"report.html"))
