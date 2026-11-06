# R-native acquisition of the observed, headerless official substrate export.
# Column adapter names follow the inspected export; downloaded text is never run.
pm_merops_source_url <- function() {
  "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/database_files/Substrate_search.txt"
}

pm_parse_merops_source <- function(source_path) {
  bytes <- readBin(source_path,"raw",n=file.info(source_path)$size)
  if (!length(bytes) || any(bytes==as.raw(0))) stop("MEROPS source is empty or contains unsupported NUL bytes.")
  # Latin-1 maps every input byte reversibly, matching the inspected Python adapter.
  text <- iconv(rawToChar(bytes),from="latin1",to="UTF-8",sub=NA)
  if (is.na(text)) stop("MEROPS source could not be decoded reversibly.")
  if (!grepl("^'CLE[0-9]+'\\t",text,perl=TRUE)) stop("MEROPS source does not begin with a supported record.")
  # Fixed splitting is linear in the source size. Repeated substring calls on
  # a long UTF-8 string can rescan preceding characters for every record.
  separator <- "\001"
  while (grepl(separator,text,fixed=TRUE)) separator <- paste0(separator,"\001")
  marked <- gsub("(?m)^('CLE[0-9]+'\\t)",paste0(separator,"\\1"),text,perl=TRUE)
  records <- strsplit(marked,separator,fixed=TRUE)[[1]]
  records <- records[nzchar(records)]
  fields <- c("merops_record_id","protease_id","substrate_name","cleavage_description",
    "p4_raw","p3_raw","p2_raw","p1_raw","p1prime_raw","p2prime_raw","p3prime_raw","p4prime_raw",
    "reference_raw","substrate_accession","p1_position","substrate_species","protease_name",
    "source_column_18_raw","source_column_19_raw","substrate_residue_range","method_code_raw",
    "additional_evidence_code_raw","cleavage_type")
  parse <- function(raw,i) {
    content <- sub("[\\r\\n]+$","",raw,perl=TRUE)
    tokens <- regmatches(content,gregexpr("(?s)('.*?'|NULL|-?[0-9]+)\\t",content,perl=TRUE))[[1]]
    if (length(tokens)!=length(fields) || !identical(paste0(tokens,collapse=""),content))
      stop("Unsupported MEROPS source structure at record ",i,"; expected 23 complete fields. No snapshot was accepted.")
    values <- substring(tokens,1L,nchar(tokens)-1L)
    quoted <- startsWith(values,"'")
    values[quoted] <- substring(values[quoted],2L,nchar(values[quoted])-1L)
    values[!quoted & values=="NULL"] <- ""
    if (!grepl("^CLE[0-9]+$",values[1]) || (nzchar(values[15]) && !grepl("^-?[0-9]+$",values[15])))
      stop("Invalid MEROPS identifier or coordinate at record ",i,".")
    values
  }
  rows <- lapply(seq_along(records),function(i) parse(records[i],i))
  if (!length(rows)) stop("MEROPS source contains no records.")
  result <- as.data.frame(do.call(rbind,rows),stringsAsFactors=FALSE)
  names(result) <- fields
  if (anyDuplicated(result$merops_record_id)) stop("MEROPS source has duplicate record IDs.")
  amino <- setNames(strsplit("ARNDCQEGHILKMFPSTWYV","")[[1]],
    strsplit("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val"," ")[[1]])
  residues <- as.matrix(result[fields[5:12]])
  canonical <- rowSums(matrix(residues %in% names(amino),nrow=nrow(result)))==8L
  encoded <- matrix(unname(amino[residues]),nrow=nrow(result),ncol=8L)
  motif <- rep("",nrow(result))
  if (any(canonical)) motif[canonical] <- apply(encoded[canonical,,drop=FALSE],1,paste0,collapse="")
  refs <- regmatches(result$reference_raw,gregexpr("\\[[^\\[\\]]+\\]",result$reference_raw,perl=TRUE))
  result$source_record_number <- seq_len(nrow(result))
  result$p4_p4prime <- motif
  result$canonical_octamer <- ifelse(canonical,"true","false")
  result$reference_ids <- vapply(refs,function(x) paste(substring(x,2L,nchar(x)-1L),collapse=";"),character(1))
  result
}

pm_fetch_merops_reference <- function(output_dir, fetch=pm_reference_fetch) {
  if (!is.function(fetch)) stop("MEROPS fetch must be a function.")
  if (!dir.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
  if (!dir.exists(output_dir) || length(list.files(output_dir,all.files=TRUE,no..=TRUE)))
    stop("MEROPS preparation requires an empty staging directory.")
  source <- file.path(output_dir,"Substrate_search.txt")
  response <- fetch(pm_merops_source_url(),source)
  if (!file.exists(source) || !is.list(response) || !identical(as.integer(response$status_code),200L))
    stop("MEROPS download did not produce a successful file response.")
  reference <- pm_parse_merops_source(source)
  snapshot <- file.path(output_dir,"merops_substrates.csv.gz")
  provenance <- file.path(output_dir,"snapshot_provenance.json")
  # Binary UTF-8 writing avoids native-codepage substitution under a C locale.
  quote_field <- function(x) paste0('"',gsub('"','""',enc2utf8(as.character(x)),fixed=TRUE),'"')
  columns <- lapply(reference,quote_field)
  lines <- c(paste(quote_field(names(reference)),collapse=","),do.call(paste,c(columns,sep=",")))
  con <- gzfile(snapshot,"wb")
  tryCatch(writeBin(charToRaw(enc2utf8(paste0(paste(lines,collapse="\n"),"\n"))),con),finally=close(con))
  headers <- response$headers
  last_modified <- if ((is.list(headers)||is.character(headers)) && "last-modified" %in% names(headers)) headers[["last-modified"]] else NULL
  retrieved <- response$retrieved_utc
  if (!is.character(retrieved) || length(retrieved)!=1L || is.na(retrieved) || !nzchar(retrieved))
    stop("MEROPS download lacks its retrieval timestamp.")
  named_counts <- function(x) {z<-table(x);as.list(setNames(as.integer(z),names(z)))}
  info <- list(schema_version="proteomaster-merops-observed-export-1",adapter_version="atelente-merops-source-1",
    source_url=pm_merops_source_url(),source_last_modified=last_modified,
    source_last_modified_index=NULL,release_version=NULL,
    version_note="The standalone substrate file does not identify its release. Source date and hashes identify this export; website release labels are not assigned to it.",
    retrieved_utc=retrieved,prepared_utc=format(Sys.time(),tz="UTC",format="%Y-%m-%dT%H:%M:%SZ"),
    source_sha256=digest::digest(file=source,algo="sha256"),source_bytes=unname(file.info(source)$size),
    source_decoding="latin-1 byte preserving; no source text is executed",
    derived_sha256=digest::digest(file=snapshot,algo="sha256"),records=nrow(reference),
    unique_record_ids=length(unique(reference$merops_record_id)),canonical_octamers=sum(reference$canonical_octamer=="true"),
    cleavage_types=named_counts(reference$cleavage_type),
    columns=lapply(seq_len(23L),function(i) list(position_1based=i,adapter_name=names(reference)[i])),
    limitations=c("Observed 23-column substrate export, not a verified SQL schema.",
      "All records and raw chemistry/evidence fields are retained; source columns 18/19 are not interpreted.",
      "Reference matching does not change candidate definitions or numerical results."))
  jsonlite::write_json(info,provenance,auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
  # A cache publisher must receive only an internally verified pair of files.
  pm_read_merops_snapshot(snapshot,provenance)
  list(snapshot=snapshot,provenance=provenance)
}
