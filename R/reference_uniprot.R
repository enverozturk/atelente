# R-native exact-accession retrieval. This layer acquires references only;
# pm_annotate_uniprot retains the full-sequence/coordinate application checks.
# Official API headers and pagination: https://doi.org/10.1093/nar/gkaf394
pm_fetch_uniprot_reference <- function(accessions, output_dir, fetch=pm_reference_fetch) {
  if (!requireNamespace("jsonlite",quietly=TRUE) || !requireNamespace("digest",quietly=TRUE))
    stop("UniProt reference retrieval requires jsonlite and digest.")
  if (!is.character(accessions) || anyNA(accessions) || any(!nzchar(accessions)))
    stop("Reference identifiers must be non-empty character values.")
  requested <- sort(unique(accessions),method="radix")
  # UniProt's documented accession syntax, retaining the existing literal
  # isoform suffix support. Do not rewrite custom/RefSeq/contaminant identifiers.
  # https://www.uniprot.org/help/accession_numbers
  queryable <- requested[grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[1-9][0-9]*)?$",requested)]
  unsupported <- setdiff(requested,queryable)
  if (!is.character(output_dir) || length(output_dir)!=1L || is.na(output_dir) || !nzchar(output_dir))
    stop("UniProt output_dir must be a directory path.")
  if (dir.exists(output_dir) && length(list.files(output_dir,all.files=TRUE,no..=TRUE)))
    stop("UniProt retrieval requires an empty staging directory.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir <- normalizePath(output_dir,winslash="/",mustWork=TRUE)
  response_dir <- file.path(output_dir,"responses")
  if (!dir.create(response_dir)) stop("Cannot create UniProt response staging directory.")
  sha <- function(path) digest::digest(file=path,algo="sha256")
  json <- function(x) as.character(jsonlite::toJSON(x,auto_unbox=TRUE,null="null",na="null",digits=I(17)))
  text_one <- function(x) is.character(x) && length(x)==1L && !is.na(x) && nzchar(x)
  write_utf8 <- function(text,path,gzip=FALSE) {
    con <- if(gzip) gzfile(path,"wb") else file(path,"wb")
    on.exit(close(con));writeBin(charToRaw(enc2utf8(text)),con)
  }
  safe_url <- function(url) text_one(url) && grepl("^https://rest\\.uniprot\\.org/uniprotkb/search\\?[^#\r\n]+$",url)
  next_page <- function(header) {
    if (is.null(header) || identical(header,"")) return(NULL)
    if (!is.character(header) || length(header)!=1L || is.na(header)) stop("Invalid UniProt pagination header.")
    parts <- strsplit(header,",(?=[[:space:]]*<)",perl=TRUE)[[1]]
    next_parts <- parts[grepl('(?:^|;)[[:space:]]*rel[[:space:]]*=[[:space:]]*(?:"next"|next)(?:[[:space:]]*;|[[:space:]]*$)',parts,perl=TRUE)]
    if (length(next_parts)>1L) stop("Multiple UniProt next-page links.")
    if (!length(next_parts)) return(NULL)
    m <- regmatches(next_parts,regexec("^[[:space:]]*<([^>]+)>",next_parts,perl=TRUE))[[1]]
    if (length(m)!=2L || !safe_url(m[2])) stop("UniProt pagination left the official HTTPS search endpoint.")
    m[2]
  }
  validate_entry <- function(entry) {
    inactive <- is.list(entry) && identical(entry$entryType,"Inactive") &&
      is.list(entry$inactiveReason) && text_one(entry$inactiveReason$inactiveReasonType)
    if (!is.list(entry) || !text_one(entry$primaryAccession) ||
        !grepl("^[A-Z0-9]+(-[0-9]+)?$",entry$primaryAccession) ||
        (!inactive && (!is.list(entry$sequence) || !text_one(entry$sequence$value) ||
        !grepl("^[A-Z]+$",entry$sequence$value))))
      stop("UniProt response entry lacks a primary accession or complete sequence.")
    if (!is.null(entry$features) && (!is.list(entry$features) ||
        any(!vapply(entry$features,is.list,logical(1)))))
      stop("UniProt response features must be an array of objects.")
    aliases <- entry$secondaryAccessions
    if (!is.null(aliases) && (!is.list(aliases) || any(!vapply(aliases,text_one,logical(1)))))
      stop("Invalid UniProt secondary accession array.")
    invisible(TRUE)
  }
  write_utf8(if(length(requested)) paste0(paste(requested,collapse="\n"),"\n") else "",file.path(output_dir,"requested_accessions.txt"))
  entries <- list();logs <- list();seen_urls <- character();release <- release_date <- NULL
  starts <- if(length(queryable)) seq.int(1L,length(queryable),by=100L) else integer()
  for (batch in seq_along(starts)) {
    ids <- queryable[starts[batch]:min(starts[batch]+99L,length(queryable))]
    query <- paste0("(",paste(paste0("accession:",ids),collapse=" OR "),")")
    url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=",utils::URLencode(query,reserved=TRUE),"&format=json&size=500")
    page <- 0L;batch_ids <- character();expected_total <- NULL
    while (!is.null(url)) {
      if (!safe_url(url) || url %in% seen_urls || page>=10000L)
        stop("Invalid or cyclic UniProt pagination.")
      seen_urls <- c(seen_urls,url)
      stem <- sprintf("batch_%03d_page_%03d",batch-1L,page)
      source <- file.path(response_dir,paste0(stem,".json"))
      meta <- fetch(url,source)
      if (!is.list(meta) || !isTRUE(meta$status_code==200L) || !file.exists(source))
        stop("UniProt download did not produce a successful response file.")
      if (!safe_url(meta$url)) stop("UniProt response URL left the official HTTPS search endpoint.")
      headers <- meta$headers
      if (is.null(names(headers))) stop("Missing UniProt HTTP headers.")
      names(headers) <- tolower(names(headers))
      if (anyDuplicated(names(headers))) stop("Ambiguous UniProt HTTP headers.")
      version <- headers[["x-uniprot-release"]];date <- headers[["x-uniprot-release-date"]]
      if (!text_one(version) || !text_one(date) || !text_one(meta$retrieved_utc))
        stop("UniProt response is missing release, release date or retrieval time.")
      if (is.null(release)) {release <- version;release_date <- date}
      else if (!identical(release,version) || !identical(release_date,date))
        stop("UniProt release changed during retrieval; refusing a mixed-release snapshot.")
      total <- headers[["x-total-results"]]
      if (!is.null(total)) {
        if (!text_one(total) || !grepl("^[0-9]+$",total)) stop("Invalid UniProt result count header.")
        total <- as.numeric(total)
        if (!is.finite(total) || (!is.null(expected_total) && total!=expected_total))
          stop("UniProt result count changed during pagination.")
        expected_total <- total
      }
      response <- tryCatch(jsonlite::fromJSON(source,simplifyVector=FALSE),error=function(e)
        stop("Invalid UniProt JSON response: ",conditionMessage(e)))
      if (!is.list(response) || !"results" %in% names(response) || !is.list(response$results) ||
          (!is.null(names(response$results)) && length(response$results)))
        stop("UniProt response requires a results array.")
      for (entry in response$results) {
        validate_entry(entry);acc <- entry$primaryAccession
        if (acc %in% batch_ids) stop("Duplicate UniProt entry within a paginated query.")
        batch_ids <- c(batch_ids,acc)
        if (!is.null(entries[[acc]]) && !identical(entries[[acc]],entry))
          stop("Inconsistent duplicate UniProt primary accession across requests: ",acc)
        entries[[acc]] <- entry
      }
      log <- list(url=url,response_url=meta$url,requested_accessions=as.list(ids),retrieved_utc=meta$retrieved_utc,
        status_code=meta$status_code,headers=headers,response_sha256=sha(source),response_bytes=unname(file.info(source)$size),
        response_file=paste0("responses/",stem,".json"))
      write_utf8(json(log),file.path(response_dir,paste0(stem,".meta.json")))
      logs[[length(logs)+1L]] <- log
      url <- next_page(headers[["link"]]);page <- page+1L
    }
    if (!is.null(expected_total) && length(batch_ids)!=expected_total)
      stop("UniProt pagination ended before the reported result count was retrieved.")
  }
  primary <- if(length(entries)) sort(names(entries),method="radix") else character()
  # Secondary accessions are audited only; the annotator still requires an exact
  # primary-accession and full-sequence match before applying any coordinates.
  alias_index <- list()
  for (acc in primary) for (alias in entries[[acc]]$secondaryAccessions)
    alias_index[[alias]] <- unique(c(alias_index[[alias]],acc))
  audit <- data.frame(requested_accession=requested,status=rep("not_returned",length(requested)),
    returned_primary_accession=rep("",length(requested)),stringsAsFactors=FALSE)
  for (i in seq_along(requested)) {
    acc <- requested[i]
    if (acc %in% unsupported) {audit$status[i] <- "unsupported_identifier";next}
    if (acc %in% primary) {
      audit$status[i] <- if(identical(entries[[acc]]$entryType,"Inactive")) "inactive" else "primary"
      audit$returned_primary_accession[i] <- acc
    }
    else if (length(alias_index[[acc]])) {
      audit$status[i] <- if(length(alias_index[[acc]])==1L) "secondary_alias" else "ambiguous_secondary_alias"
      audit$returned_primary_accession[i] <- paste(sort(alias_index[[acc]],method="radix"),collapse=";")
    }
  }
  # Custom identifiers may be Unicode even in a C-locale R session.
  csv_field <- function(x) if(!length(x)) character() else paste0('"',gsub('"','""',enc2utf8(x),fixed=TRUE),'"')
  audit_lines <- c(paste(csv_field(names(audit)),collapse=","),
    do.call(paste,c(lapply(audit,csv_field),sep=",")))
  write_utf8(paste0(paste(audit_lines,collapse="\n"),"\n"),file.path(output_dir,"request_audit.csv"))
  payload <- json(list(results=unname(entries[primary])))
  snapshot <- file.path(output_dir,"uniprot_processing_records.json.gz")
  write_utf8(payload,snapshot,gzip=TRUE)
  retrieved <- vapply(logs,function(x) x$retrieved_utc,character(1))
  provenance <- list(schema_version="proteomaster-uniprot-rest-1",derived_sha256=sha(snapshot),records=length(entries),
    uncompressed_sha256=digest::digest(charToRaw(enc2utf8(payload)),algo="sha256",serialize=FALSE),
    uncompressed_bytes=length(charToRaw(enc2utf8(payload))),compressed_bytes=unname(file.info(snapshot)$size),
    release=release,release_date=release_date,retrieved_utc_start=if(length(retrieved)) min(retrieved) else NULL,
    retrieved_utc_end=if(length(retrieved)) max(retrieved) else NULL,prepared_utc=format(Sys.time(),tz="UTC",format="%Y-%m-%dT%H:%M:%SZ"),
    requested_accessions=length(requested),requested_accession_ids=as.list(requested),
    queried_accession_ids=as.list(queryable),unsupported_accession_ids=as.list(unsupported),
    exact_primary_matches=sum(audit$status=="primary"),secondary_alias_matches=sum(audit$status=="secondary_alias"),
    ambiguous_secondary_aliases=sum(audit$status=="ambiguous_secondary_alias"),
    inactive_accessions=as.list(audit$requested_accession[audit$status=="inactive"]),
    not_returned=as.list(audit$requested_accession[audit$status=="not_returned"]),
    request_audit_sha256=sha(file.path(output_dir,"request_audit.csv")),
    scope="Exact accession queries from mapped records; full JSON entries retained. Secondary aliases and missing accessions are audited, never automatically transferred. No peptide sequences or quantities are sent.",
    retrieval_status=if(!length(requested)) "empty_request" else if(!length(queryable)) "no_queryable_accessions" else if(length(unsupported)) "partial_identifier_coverage" else "complete",requests=logs)
  prov_path <- file.path(output_dir,"snapshot_provenance.json")
  write_utf8(json(provenance),prov_path)
  # Validate the same on-disk contract the annotation/cache layer consumes.
  if (exists("pm_read_uniprot_snapshot",mode="function")) pm_read_uniprot_snapshot(snapshot,prov_path)
  list(snapshot=snapshot,provenance=prov_path)
}
