# Repackage an existing report without recomputing mapping, quantities or tests.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: Rscript tools/refresh_output_bundle.R INPUT_REPORT NEW_DIRECTORY")
if(.Platform$OS.type=="windows")invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8")))
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))source(f,encoding="UTF-8")
refresh_bundle <- function(input_report, destination) {
input <- normalizePath(input_report,winslash="/",mustWork=TRUE)
out <- destination
if(file.exists(out))stop("Choose a new output directory.")
source_root <- dirname(input)
html <- paste(readLines(input,encoding="UTF-8",warn=FALSE),collapse="\n")
marker <- '<script type="application/json" id="data">'
at <- regexpr(marker,html,fixed=TRUE)[1]
if(at<1L)stop("Report data block not found.")
tail <- substring(html,at+nchar(marker),nchar(html));end<-regexpr('</script>',tail,fixed=TRUE)[1]
if(end<1L)stop("Report data block is incomplete.")
# Match the existing single-file transport ceiling. This prevents the known
# oversized legacy conversion from exhausting RAM; it is not a RAM guarantee
# for every smaller report. A chunked legacy importer is a separate workflow.
payload <- pm_report_decode_json(substr(tail,1,end-1L),max_uncompressed_bytes=450*1024^2)
original <- payload
manifest_path <- file.path(source_root,"run_manifest.json")
manifest <- if(file.exists(manifest_path))jsonlite::read_json(manifest_path,simplifyVector=FALSE) else list(tool="Atelente")
hash <- function(p)digest::digest(file=p,algo="sha256")
# Older wizard reports kept source-column/scale metadata only in their audit.
# Restore it only from the exact generated artifact, with a recorded hash and
# matching report columns. Never infer source headers from canonical aliases.
if (is.null(payload$column_mapping_audit)) {
  relative <- c("audit/column_mapping_audit.json", "column_mapping_audit.json")
  available <- relative[file.exists(file.path(source_root,relative))]
  if (length(available)>1L) stop("Multiple column mapping audits found; choose an unambiguous source bundle.")
  if (length(available)) {
    relative <- available[1]; mapping_path <- file.path(source_root,relative)
    if (!startsWith(normalizePath(mapping_path,winslash="/",mustWork=TRUE),paste0(source_root,"/")))
      stop("Column mapping audit escapes the source report directory.")
    expected <- manifest$output_sha256[[relative]]
    if (identical(relative,"audit/column_mapping_audit.json")) {
      index_path <- file.path(source_root,"audit","audit_files.json")
      if (file.exists(index_path)) {
        if (!is.null(manifest$audit_index_sha256) && !identical(hash(index_path),manifest$audit_index_sha256))
          stop("Column mapping audit index SHA256 mismatch.")
        audit_index <- jsonlite::read_json(index_path,simplifyVector=FALSE)
        if (!identical(audit_index$schema,"atelente_audit_v1")) stop("Invalid column mapping audit index schema.")
        indexed <- audit_index$sha256[["column_mapping_audit.json"]]
        if (!is.null(indexed)) {
          if (!is.null(expected) && !identical(expected,indexed)) stop("Conflicting column mapping audit hashes.")
          expected <- indexed
        }
      }
    }
    if (!is.character(expected)||length(expected)!=1L||is.na(expected)||!grepl("^[a-fA-F0-9]{64}$",expected))
      stop("Column mapping audit has no verified SHA256 record; restore its original manifest or rerun the analysis.")
    if (!identical(hash(mapping_path),tolower(expected))) stop("Column mapping audit SHA256 mismatch.")
    mapping <- jsonlite::read_json(mapping_path,simplifyVector=TRUE)
    q <- mapping$quantities; samples <- payload$pair_view$samples
    columns <- if(is.data.frame(samples)&&"Peptide column" %in% names(samples))samples[["Peptide column"]] else character()
    # JSON [] from a mapping-only run has no data-frame column attributes.
    if(is.list(q)&&!length(q)&&!length(columns))
      q<-data.frame(source_column=character(),canonical_column=character(),input_scale=character())
    if (!identical(mapping$schema,"atelente-column-mapping-1") || !is.data.frame(q) ||
        !all(c("source_column","canonical_column","input_scale") %in% names(q)) ||
        !all(vapply(q[c("source_column","canonical_column","input_scale")],is.character,logical(1))) ||
        anyNA(q[c("source_column","canonical_column","input_scale")]) ||
        any(!nzchar(q$source_column)) || any(!nzchar(q$canonical_column)) ||
        anyDuplicated(q$source_column) || anyDuplicated(q$canonical_column) ||
        any(!q$input_scale %in% c("linear","log2")))
      stop("Invalid column mapping audit metadata.")
    if(anyNA(columns)||anyDuplicated(columns)||!setequal(q$canonical_column,columns))
      stop("Column mapping audit canonical columns do not match report samples.")
    payload$column_mapping_audit <- mapping
  }
}
manifest$presentation_refresh <- list(source_report_sha256=hash(input),version="reader_export_v7",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
# Only generated artifacts in the source report directory; no input/reference fetch.
files <- list.files(source_root,full.names=TRUE,recursive=TRUE,all.files=TRUE,no..=TRUE)
files <- files[!dir.exists(files)&grepl("\\.(csv|json|xlsx)$",files,ignore.case=TRUE)]
# Rebuild presentation artifacts in the new directory, retaining original analyses.
files <- files[!basename(files) %in% c("atelente_results.xlsx","atelente_tables.xlsx","audit_files.json")]
inside <- function(p)startsWith(normalizePath(p,winslash="/",mustWork=TRUE),paste0(source_root,"/"))
if(any(!vapply(files,inside,logical(1))))stop("Generated artifact escapes source directory.")
# Build beside the requested destination and publish only a complete bundle.
# on.exit also removes the staged copy after export/identity failures.
parent <- dirname(out)
dir.create(parent,recursive=TRUE,showWarnings=FALSE)
parent <- normalizePath(parent,winslash="/",mustWork=TRUE)
destination <- file.path(parent,basename(out))
out <- tempfile(".atelente-repack-",tmpdir=parent)
if(!dir.create(out))stop("Could not create temporary output directory.")
staging <- normalizePath(out,winslash="/",mustWork=TRUE)
if(!startsWith(staging,paste0(parent,"/")))stop("Temporary output escapes its parent directory.")
on.exit(if(dir.exists(staging))unlink(staging,recursive=TRUE),add=TRUE)
for(src in files){relative<-substring(gsub("\\\\","/",src),nchar(source_root)+2L);dest<-file.path(out,relative);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);if(!file.copy(src,dest)||hash(src)!=hash(dest))stop("Artifact copy failed.")}
payload <- pm_prepare_export_bundle(payload,out,manifest,TRUE)
keep <- setdiff(names(original),c("excel_file","sample_metadata_file","analysis_files","export_workbook","audit_file","export_settings","pair_view","quantity_files","evidence_workbook"))
stopifnot(identical(payload[keep],original[keep]))
prefix <- function(path) if(startsWith(path,"audit/"))path else paste0("audit/",path)
if(length(original$quantity_files))stopifnot(identical(payload$quantity_files,lapply(original$quantity_files,prefix)))
if(!is.null(original$evidence_workbook))stopifnot(identical(payload$evidence_workbook,prefix(original$evidence_workbook)))
pv <- payload$pair_view;ov <- original$pair_view;pv$selection<-ov$selection<-NULL;pv$selection_k<-ov$selection_k<-NULL
# Legacy reports gain derived R summaries on their first compact refresh. An
# existing prepared summary must remain unchanged on subsequent refreshes.
if(is.null(ov$prepared_summaries)) pv$prepared_summaries <- NULL
if(is.null(ov$spanning_end_reporting)) pv$spanning_end_reporting <- NULL
stopifnot(identical(pv,ov))
pm_write_report(payload,file.path(root,"inst","templates","report_template.html"),file.path(out,"report.html"))
manifest <- pm_finalize_audit(out,manifest)
for(path in unlist(payload$quantity_files,use.names=FALSE)) {
  actual <- normalizePath(file.path(out,path),winslash="/",mustWork=TRUE)
  if(!startsWith(actual,paste0(staging,"/")))stop("Quantity download escapes the output directory.")
}
manifest$output_sha256 <- as.list(setNames(vapply(file.path(out,c("report.html","atelente_results.xlsx","atelente_audit.zip")),hash,character(1)),c("report.html","atelente_results.xlsx","atelente_audit.zip")))
jsonlite::write_json(manifest,file.path(out,"run_manifest.json"),auto_unbox=TRUE,pretty=TRUE,na="null",null="null",digits=I(17))
if(file.exists(destination)||!file.rename(out,destination))stop("Could not publish completed output directory.")
out <- destination
cat(jsonlite::toJSON(list(directory=normalizePath(out,winslash="/"),report_bytes=file.info(file.path(out,"report.html"))$size,workbook_bytes=file.info(file.path(out,"atelente_results.xlsx"))$size,audit_bytes=file.info(file.path(out,"atelente_audit.zip"))$size),auto_unbox=TRUE),"\n")
}
refresh_bundle(args[1],args[2])
