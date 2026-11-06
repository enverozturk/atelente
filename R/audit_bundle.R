# External audit packaging; no analysis or reference downloads occur here.
pm_finalize_audit <- function(output_dir, manifest) {
  for(package in c("digest","jsonlite","zip"))
    if(!requireNamespace(package,quietly=TRUE)) stop("Install ",package," for audit packaging.")
  root <- normalizePath(output_dir,winslash="/",mustWork=TRUE)
  audit <- file.path(root,"audit")
  dir.create(audit,showWarnings=FALSE)
  inside <- function(path) startsWith(paste0(normalizePath(path,winslash="/",mustWork=TRUE),"/"),paste0(root,"/"))
  if(!inside(audit)) stop("Audit directory escapes output directory.")
  hash <- function(path) digest::digest(file=path,algo="sha256")
  # Keep the original manifest beside its original artifacts. The root manifest
  # can then describe the new package layout without destroying input hashes.
  original <- file.path(root,"run_manifest.json")
  audit_manifest <- file.path(audit,"run_manifest.json")
  if(!file.exists(audit_manifest)) {
    if(file.exists(original)) {
      if(!file.copy(original,audit_manifest,overwrite=FALSE)||hash(original)!=hash(audit_manifest)) stop("Audit manifest copy failed.")
    } else jsonlite::write_json(manifest,audit_manifest,pretty=TRUE,auto_unbox=TRUE,na="null",null="null",digits=I(17))
  }
  paths <- list.files(root,pattern="\\.(csv|json|xlsx)$",full.names=TRUE,ignore.case=TRUE)
  paths <- paths[!basename(paths) %in% c("atelente_results.xlsx","run_manifest.json")]
  moved <- list()
  # Historical quantity-attachment reports keep generated tables below these
  # directories. Include them in the full audit without collecting arbitrary
  # input/reference directories that may sit beside a report.
  for(name in c("quantities","statistics")) {
    source <- file.path(root,name)
    if(!dir.exists(source)) next
    target <- file.path(audit,name)
    if(file.exists(target)) stop("Audit artifact already exists: ",name)
    if(nzchar(Sys.readlink(source))||!inside(source)) stop("Unsafe generated audit directory: ",name)
    # Check directory links before descending, preventing junction traversal.
    pending <- source; contents <- character()
    while(length(pending)) {
      directory <- pending[1];pending <- pending[-1]
      children <- list.files(directory,full.names=TRUE,all.files=TRUE,no..=TRUE)
      for(child in children) {
        if(nzchar(Sys.readlink(child))||!inside(child)) stop("Unsafe generated audit entry: ",name)
        if(dir.exists(child))pending<-c(pending,child) else contents<-c(contents,child)
      }
    }
    relative <- substring(gsub("\\\\","/",contents),nchar(source)+2L)
    expected <- vapply(contents,hash,character(1))
    if(!file.rename(source,target)) stop("Audit directory move failed: ",name)
    actual <- vapply(file.path(target,relative),hash,character(1))
    if(!identical(unname(expected),unname(actual))) stop("Audit directory hash mismatch: ",name)
    moved[[name]] <- paste0("audit/",name)
  }
  for(source in paths) {
    if(!inside(source)||dir.exists(source)) stop("Unsafe audit artifact.")
    target <- file.path(audit,basename(source))
    if(file.exists(target)) stop("Audit artifact already exists: ",basename(source))
    expected <- hash(source)
    if(!file.rename(source,target)||hash(target)!=expected) stop("Audit artifact move failed: ",basename(source))
    moved[[basename(source)]] <- paste0("audit/",basename(source))
  }
  readme <- c("Atelente full audit", "", "These are generated analysis records, not a second analysis.",
    "Original CSV values, IDs, source-row links, exclusions, settings and reference provenance are retained.",
    "The original run_manifest.json describes the original artifact names and hashes.",
    "audit_files.json hashes the files in this audit directory; paths are relative to audit/.",
    "The root run_manifest.json describes the packaged outputs and ZIP hash.",
    "External FASTA/reference databases and original input workbooks are not copied here.",
    "Re-running from raw inputs also requires those inputs/references, identified by the recorded hashes.",
    "Keep this audit directory with report.html when refreshing a report.")
  writeLines(readme,file.path(audit,"README.txt"),useBytes=TRUE)
  files <- list.files(audit,recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE)
  files <- files[!dir.exists(files)&basename(files)!="audit_files.json"]
  if(any(!vapply(files,inside,logical(1)))) stop("Audit file escapes output directory.")
  relative <- substring(gsub("\\\\","/",files),nchar(audit)+2L)
  hashes <- stats::setNames(lapply(files,hash),relative)
  jsonlite::write_json(list(schema="atelente_audit_v1",sha256=hashes),file.path(audit,"audit_files.json"),pretty=TRUE,auto_unbox=TRUE)
  archive <- file.path(root,"atelente_audit.zip")
  zip::zipr(archive,files="audit",root=root,include_directories=TRUE)
  if(!file.exists(archive)) stop("Audit archive was not created.")
  manifest$artifact_layout <- "atelente_bundle_v1"
  manifest$audit_file <- "atelente_audit.zip"
  manifest$audit_sha256 <- hash(archive)
  manifest$audit_paths <- moved
  manifest$audit_index_sha256 <- hash(file.path(audit,"audit_files.json"))
  manifest
}
