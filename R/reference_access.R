# R-native reference transport and immutable, validated local cache.
pm_reference_fetch <- function(url, path) {
  if (!is.character(url) || length(url)!=1L || !grepl("^https://(ftp\\.ebi\\.ac\\.uk|rest\\.uniprot\\.org)/",url))
    stop("Reference download requires an official HTTPS endpoint.")
  if (!requireNamespace("curl",quietly=TRUE)) stop("Install the curl R package for reference downloads.")
  last <- NULL
  for (attempt in seq_len(3L)) {
    response <- tryCatch(curl::curl_fetch_disk(url,path,handle=curl::new_handle(
      timeout=120,connecttimeout=30,followlocation=TRUE,maxredirs=5,
      useragent="Atelente-reference-client/1",accept_encoding="")),error=function(e)e)
    if (!inherits(response,"error") && identical(as.integer(response$status_code),200L))
      return(list(headers=curl::parse_headers_list(response$headers),status_code=200L,
        url=response$url,retrieved_utc=format(Sys.time(),tz="UTC",format="%Y-%m-%dT%H:%M:%SZ")))
    last <- if (inherits(response,"error")) conditionMessage(response) else paste("HTTP",response$status_code)
    if (!inherits(response,"error") && !response$status_code %in% c(408,429,500,502,503,504)) break
    if (attempt<3L) Sys.sleep(2^(attempt-1L))
  }
  stop("Reference download failed: ",last,". Retry with connectivity, supply a verified snapshot, or explicitly set annotations = FALSE.")
}

pm_reference_cache_dir <- function(cache_dir=NULL) {
  if (is.null(cache_dir)) cache_dir <- tools::R_user_dir("atelente","cache")
  if (!is.character(cache_dir)||length(cache_dir)!=1L||is.na(cache_dir)||!nzchar(cache_dir))
    stop("cache_dir must be a non-empty directory path.")
  dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  if (!dir.exists(cache_dir)) stop("Cannot create reference cache directory: ",cache_dir)
  normalizePath(cache_dir,winslash="/",mustWork=TRUE)
}

pm_reference_validate <- function(provider,snapshot,provenance) {
  loaded <- switch(provider,merops=pm_read_merops_snapshot(snapshot,provenance),
    uniprot=pm_read_uniprot_snapshot(snapshot,provenance),stop("Unknown reference provider."))
  loaded$provenance
}

pm_reference_cache_read <- function(directory,provider,requested=NULL) {
  meta_path <- file.path(directory,"cache_meta.json")
  if (!file.exists(meta_path)) stop("Incomplete reference cache: ",directory,". Supply a snapshot or use a new cache directory.")
  meta <- jsonlite::read_json(meta_path,simplifyVector=TRUE)
  if (!identical(meta$schema,"atelente-reference-cache-1")||!identical(meta$provider,provider))
    stop("Invalid reference cache metadata: ",directory)
  # Only literal basenames generated here may be resolved inside a cache entry.
  safe <- function(name) is.character(name)&&length(name)==1L&&!is.na(name)&&
    grepl("^[A-Za-z0-9_.-]+$",name)&&!name %in% c(".","..")
  if (!safe(meta$snapshot)||!safe(meta$provenance)) stop("Invalid cached reference filenames.")
  snapshot <- file.path(directory,meta$snapshot);provenance <- file.path(directory,meta$provenance)
  sha <- function(path) digest::digest(file=path,algo="sha256")
  if (!file.exists(snapshot)||!file.exists(provenance)||
      !identical(sha(snapshot),meta$snapshot_sha256)||!identical(sha(provenance),meta$provenance_sha256))
    stop("Reference cache integrity check failed: ",directory,". Supply a verified snapshot or use a new cache directory.")
  ids <- as.character(unlist(meta$requested_accessions,use.names=FALSE))
  if (!is.null(requested)&&!all(requested %in% ids)) stop("Reference cache does not cover requested accessions.")
  prov <- pm_reference_validate(provider,snapshot,provenance)
  if (provider=="uniprot"&&!setequal(ids,as.character(unlist(prov$requested_accession_ids,use.names=FALSE))))
    stop("Cached request coverage disagrees with the verified UniProt provenance.")
  list(snapshot=snapshot,provenance=provenance,metadata=prov,mode="cache")
}

pm_reference_resolve <- function(provider,spec=NULL,cache_dir=NULL,accessions=character(),fetch=pm_reference_fetch) {
  if (!is.null(spec)) {
    if (!is.list(spec)||!is.character(spec$snapshot)||length(spec$snapshot)!=1L||
        !is.character(spec$provenance)||length(spec$provenance)!=1L)
      stop(provider," requires snapshot and provenance paths.")
    snapshot <- normalizePath(spec$snapshot,winslash="/",mustWork=TRUE)
    provenance <- normalizePath(spec$provenance,winslash="/",mustWork=TRUE)
    return(list(snapshot=snapshot,provenance=provenance,
      metadata=pm_reference_validate(provider,snapshot,provenance),mode="supplied snapshot"))
  }
  cache <- pm_reference_cache_dir(cache_dir)
  ids <- unique(as.character(accessions))
  if (anyNA(ids)||any(!nzchar(ids))) stop("Reference accessions must be non-empty.")
  ids <- sort(ids,method="radix")
  key <- if (provider=="merops") "substrates-v1" else paste0("accessions-v1-",digest::digest(paste(ids,collapse="\n"),algo="sha256",serialize=FALSE))
  base <- file.path(cache,provider);dir.create(base,recursive=TRUE,showWarnings=FALSE)
  target <- file.path(base,key)
  if (dir.exists(target)) return(pm_reference_cache_read(target,provider,if(provider=="uniprot") ids else NULL))
  # A cached superset also answers a subset offline, including recorded no-return IDs.
  if (provider=="uniprot") {
    dirs <- sort(list.dirs(base,recursive=FALSE,full.names=TRUE))
    dirs <- dirs[startsWith(basename(dirs),"accessions-v1-")]
    for (dir in dirs) {
      meta_file <- file.path(dir,"cache_meta.json")
      if (!file.exists(meta_file)) next
      meta <- tryCatch(jsonlite::read_json(meta_file,simplifyVector=TRUE),error=function(e)NULL)
      covered <- as.character(unlist(meta$requested_accessions,use.names=FALSE))
      if (all(ids %in% covered)) return(pm_reference_cache_read(dir,provider,ids))
    }
  }
  stage <- tempfile(pattern=".stage-",tmpdir=base)
  dir.create(stage)
  # No partial directory is published. Failed stages are diagnostic remnants,
  # ignored by future lookups; never delete an unrelated user cache directory.
  message("Downloading ",provider," reference (first use for this request)...")
  built <- switch(provider,merops=pm_fetch_merops_reference(stage,fetch=fetch),
    uniprot=pm_fetch_uniprot_reference(ids,stage,fetch=fetch),stop("Unknown reference provider."))
  stage_real <- normalizePath(stage,winslash="/",mustWork=TRUE)
  for (field in c("snapshot","provenance")) {
    actual <- normalizePath(built[[field]],winslash="/",mustWork=TRUE)
    if (!identical(dirname(actual),stage_real)) stop("Reference provider wrote outside its staging directory.")
  }
  prov <- pm_reference_validate(provider,built$snapshot,built$provenance)
  sha <- function(path) digest::digest(file=path,algo="sha256")
  meta <- list(schema="atelente-reference-cache-1",provider=provider,
    snapshot=basename(built$snapshot),provenance=basename(built$provenance),
    snapshot_sha256=sha(built$snapshot),provenance_sha256=sha(built$provenance),
    requested_accessions=if(provider=="uniprot") ids else character(),
    created_utc=format(Sys.time(),tz="UTC",format="%Y-%m-%dT%H:%M:%SZ"))
  jsonlite::write_json(meta,file.path(stage,"cache_meta.json"),auto_unbox=TRUE,pretty=TRUE)
  # Concurrent first-use writers may race; a verified winner can be reused.
  if (!dir.exists(target)&&file.rename(stage,target)) {
    resolved <- pm_reference_cache_read(target,provider,if(provider=="uniprot") ids else NULL)
    resolved$mode <- if(identical(prov$retrieval_status,"no_queryable_accessions")) "local identifier audit" else "downloaded";return(resolved)
  }
  if (dir.exists(target)) return(pm_reference_cache_read(target,provider,if(provider=="uniprot") ids else NULL))
  stop("Cannot publish the validated reference cache: ",target)
}

pm_reference_info <- function(resolved,provider) {
  p <- resolved$metadata
  first <- function(...) {for (x in list(...)) if(!is.null(x)&&length(x)&&!is.na(x[1])&&nzchar(as.character(x[1]))) return(as.character(x[1]));NULL}
  list(provider=provider,access=resolved$mode,
    release=first(p$release,p$release_version),
    source_date=first(p$release_date,p$source_last_modified,p$source_last_modified_index),
    retrieved_utc=first(p$retrieved_utc,p$retrieved_utc_end,p$prepared_utc),
    snapshot_sha256=digest::digest(file=resolved$snapshot,algo="sha256"),
    provenance_sha256=digest::digest(file=resolved$provenance,algo="sha256"),
    records=p$records,source_url=p$source_url,
    unsupported_accession_ids=p$unsupported_accession_ids)
}

pm_annotate_references <- function(result,annotations=TRUE,merops=NULL,uniprot=NULL,cache_dir=NULL) {
  if(!is.logical(annotations)||length(annotations)!=1L||is.na(annotations)) stop("annotations must be TRUE or FALSE.")
  if(!annotations){
    if(!is.null(merops)||!is.null(uniprot)) stop("annotations=FALSE conflicts with supplied reference snapshots.")
    result$merops <- result$uniprot <- NULL
    result$reference_access <- list(enabled=FALSE,status="disabled",message="Annotations disabled by user (annotations = FALSE).")
    return(result)
  }
  providers <- list()
  if(!is.null(result$cleavage)) {
    r <- pm_reference_resolve("merops",merops,cache_dir)
    result <- pm_annotate_merops(result,r$snapshot,r$provenance)
    providers$merops <- pm_reference_info(r,"merops")
  } else if(!is.null(merops)) stop("MEROPS site annotations require a cleavage profile.")
  r <- pm_reference_resolve("uniprot",uniprot,cache_dir,accessions=unique(result$matches$accession))
  result <- pm_annotate_uniprot(result,r$snapshot,r$provenance)
  providers$uniprot <- pm_reference_info(r,"uniprot")
  unsupported <- as.character(unlist(r$metadata$unsupported_accession_ids,use.names=FALSE))
  unsupported <- intersect(unsupported,unique(result$matches$accession))
  result$reference_access <- list(enabled=TRUE,status=if(length(unsupported)) "partial_identifier_coverage" else "complete",providers=providers,
    unsupported_uniprot_ids=as.list(unsupported),
    merops_applicability=if(is.null(result$cleavage)) "No cleavage layer; site annotation not applicable" else "Site annotation completed")
  result
}
