# Optional descriptive processing annotation; never filters or reclassifies observations.
pm_up_text <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "" else as.character(x)
pm_up_json <- function(x) as.character(jsonlite::toJSON(x, auto_unbox=TRUE, null="null", na="null"))
pm_up_rows <- function(rows, empty) if (length(rows)) { x <- do.call(rbind,rows);rownames(x)<-NULL;x } else empty

pm_read_uniprot_snapshot <- function(snapshot_path, provenance_path) {
  if (!requireNamespace("digest",quietly=TRUE) || !requireNamespace("jsonlite",quietly=TRUE))
    stop("Install digest and jsonlite for UniProt annotation.")
  provenance <- jsonlite::fromJSON(provenance_path,simplifyVector=FALSE)
  if (!identical(provenance$schema_version,"proteomaster-uniprot-rest-1")) stop("Unsupported UniProt snapshot schema_version.")
  hash <- provenance$derived_sha256
  if (!is.character(hash) || length(hash)!=1L || !grepl("^[0-9a-fA-F]{64}$",hash)) stop("UniProt provenance requires derived_sha256.")
  actual <- digest::digest(file=snapshot_path,algo="sha256")
  if (!identical(tolower(hash),actual)) stop("UniProt snapshot SHA-256 does not match provenance.")
  # UTF-8 text connections may convert to the native C locale and truncate.
  # Read the exact decompressed bytes, then mark JSON text explicitly as UTF-8.
  con <- if(grepl("\\.gz$",snapshot_path,ignore.case=TRUE)) gzfile(snapshot_path,"rb") else file(snapshot_path,"rb")
  on.exit(close(con))
  chunks <- list()
  repeat {
    bytes <- readBin(con,what="raw",n=1048576L)
    if(!length(bytes)) break
    chunks[[length(chunks)+1L]] <- bytes
  }
  text <- rawToChar(if(length(chunks)) do.call(c,chunks) else raw())
  Encoding(text) <- "UTF-8"
  snapshot <- jsonlite::fromJSON(text,simplifyVector=FALSE)
  entries <- snapshot$results
  if (!is.list(entries) || is.null(provenance$records) || length(provenance$records)!=1L ||
      !is.numeric(provenance$records) || provenance$records!=length(entries)) stop("UniProt snapshot requires a results array and matching records count.")
  if (any(!vapply(entries,function(e) is.list(e) && nzchar(pm_up_text(e$primaryAccession)),logical(1))))
    stop("Each UniProt result requires a primaryAccession.")
  provenance$verified_snapshot_sha256 <- actual
  provenance$verification <- "Snapshot bytes and entry count verified against supplied provenance; provenance is not digitally authenticated."
  list(entries=entries,provenance=provenance)
}

pm_up_position <- function(position) {
  value <- if (is.list(position)) position$value else NULL
  modifier <- if (is.list(position)) pm_up_text(position$modifier) else ""
  valid <- is.numeric(value) && length(value)==1L && !is.na(value) && is.finite(value) && value==floor(value) && value>=1 && value<=.Machine$integer.max
  list(value=if(valid) as.integer(value) else NA_integer_,modifier=modifier,exact=valid && identical(modifier,"EXACT"))
}

pm_annotate_uniprot <- function(result,snapshot_path,provenance_path) {
  loaded <- pm_read_uniprot_snapshot(snapshot_path,provenance_path)
  entries <- loaded$entries
  f <- result$fasta;m <- result$matches
  if (!is.data.frame(f) || !all(c("record_id","accession","sequence") %in% names(f)) ||
      !is.data.frame(m) || !all(c("record_id","peptide_id","start","end") %in% names(m))) stop("UniProt annotation requires FASTA records and peptide matches.")
  primary <- vapply(entries,function(e) pm_up_text(e$primaryAccession),character(1))
  unsupported <- as.character(unlist(loaded$provenance$unsupported_accession_ids,use.names=FALSE))
  by_accession <- split(seq_along(entries),primary)
  relevant <- f[f$record_id %in% unique(m$record_id),,drop=FALSE]
  supported <- c("Signal","Transit peptide","Propeptide","Chain","Peptide","Initiator methionine")
  record_rows <- feature_rows <- boundary_rows <- list();ri<-fi<-bi<-0L
  for (r in seq_len(nrow(relevant))) {
    fr <- relevant[r,,drop=FALSE];accession<-as.character(fr$accession);sequence<-as.character(fr$sequence)
    idx <- by_accession[[accession]]
    status <- if (!length(idx)) {if(accession %in% unsupported) "unsupported_identifier" else if(grepl("-[0-9]+$",accession)) "unsupported_isoform" else "not_found"} else
      if(length(idx)>1L) "ambiguous_accession" else if(!nzchar(pm_up_text(entries[[idx]]$sequence$value))) "sequence_unavailable" else
        if(!identical(sequence,pm_up_text(entries[[idx]]$sequence$value))) "sequence_mismatch" else "matched"
    one <- if(length(idx)==1L) entries[[idx]] else list()
    ri<-ri+1L
    record_rows[[ri]] <- data.frame(record_id=fr$record_id,accession=accession,status=status,
      uniprot_accession=if(length(idx)) paste(primary[idx],collapse=";") else "",
      entry_version=pm_up_text(one$entryAudit$entryVersion),sequence_version=pm_up_text(one$entryAudit$sequenceVersion),
      annotation_date=pm_up_text(one$entryAudit$lastAnnotationUpdateDate),sequence_date=pm_up_text(one$entryAudit$lastSequenceUpdateDate),
      fasta_sequence_sha256=digest::digest(sequence,algo="sha256",serialize=FALSE),
      uniprot_sequence_sha256=if(nzchar(pm_up_text(one$sequence$value))) digest::digest(one$sequence$value,algo="sha256",serialize=FALSE) else "",
      stringsAsFactors=FALSE)
    for (ei in idx) {
      entry<-entries[[ei]]
      for(j in seq_along(entry$features)) {
        ft<-entry$features[[j]];type<-pm_up_text(ft$type)
        if(!type %in% supported) next
        a<-pm_up_position(ft$location$start);b<-pm_up_position(ft$location$end)
        sequence_context<-pm_up_text(ft$location$sequence)
        # UniProt location.sequence identifies a distinct coordinate sequence.
        context_ok<-!nzchar(sequence_context) || identical(sequence_context,accession)
        ordered<-is.na(a$value) || is.na(b$value) || a$value<=b$value
        ae<-a$exact && a$value<=nchar(sequence);be<-b$exact && b$value<=nchar(sequence)
        compatible<-status=="matched" && context_ok && ordered
        cs<-if(status!="matched") status else if(!context_ok) "different_sequence_context" else if(!ordered) "invalid_range" else
          if((a$exact && !ae)||(b$exact && !be)) "out_of_range" else if(ae&&be) "compatible" else if(ae||be) "partial_exact_boundaries" else "uncertain_boundaries"
        if(cs=="out_of_range") compatible<-FALSE
        description<-pm_up_text(ft$description)
        evidences<-if(is.null(ft$evidences)) list() else ft$evidences
        codes<-unique(vapply(evidences,function(e) pm_up_text(e$evidenceCode),character(1)))
        evidence_codes<-paste(codes[nzchar(codes)],collapse=";")
        key<-paste(fr$record_id,ei,j,sep=":")
        fi<-fi+1L
        feature_rows[[fi]]<-data.frame(feature_key=key,record_id=fr$record_id,accession=accession,
          uniprot_accession=primary[ei],entry_index=ei,feature_index=j,feature_type=type,description=description,
          source_feature_id=pm_up_text(ft$featureId),start=a$value,end=b$value,start_modifier=a$modifier,end_modifier=b$modifier,
          sequence_context=sequence_context,coordinate_status=cs,evidence_codes=evidence_codes,evidence_json=pm_up_json(evidences),
          location_json=pm_up_json(ft$location),stringsAsFactors=FALSE)
        if(!compatible) next
        ks<-integer();labels<-character()
        if(type %in% c("Chain","Peptide")) {
          if(ae){ks<-c(ks,a$value-1L);labels<-c(labels,paste0(tolower(type),"_start"))}
          if(be){ks<-c(ks,b$value);labels<-c(labels,paste0(tolower(type),"_end"))}
        } else if(type=="Propeptide") {
          if(ae && a$value>1L){ks<-c(ks,a$value-1L);labels<-c(labels,"before_propeptide")}
          if(be && b$value<nchar(sequence)){ks<-c(ks,b$value);labels<-c(labels,"after_propeptide")}
        } else if(type %in% c("Signal","Transit peptide")) {
          if(be && b$value<nchar(sequence)){ks<-b$value;labels<-if(type=="Signal") "after_signal" else "after_transit_peptide"}
        } else if(type=="Initiator methionine") {
          # Do not turn retained or unrecognized initiation annotations into cleavage.
          if(ae && be && a$value==1L && b$value==1L && startsWith(sequence,"M") && tolower(trimws(description)) %in% c("removed","removed; alternate") && nchar(sequence)>1L) {
            ks<-1L;labels<-"initiator_methionine_removed"
          }
        }
        for(k in seq_along(ks)) {
          bi<-bi+1L
          boundary_rows[[bi]]<-data.frame(record_id=fr$record_id,accession=accession,bond_left=ks[k],bond_right=ks[k]+1L,
            feature_key=key,feature_type=type,description=description,boundary=labels[k],
            internal_bond=ks[k]>0L && ks[k]<nchar(sequence),evidence_codes=evidence_codes,evidence_json=pm_up_json(evidences),
            uniprot_accession=primary[ei],source_feature_id=pm_up_text(ft$featureId),stringsAsFactors=FALSE)
        }
      }
    }
  }
  er<-data.frame(record_id=integer(),accession=character(),status=character(),uniprot_accession=character(),entry_version=character(),sequence_version=character(),annotation_date=character(),sequence_date=character(),fasta_sequence_sha256=character(),uniprot_sequence_sha256=character())
  ef<-data.frame(feature_key=character(),record_id=integer(),accession=character(),uniprot_accession=character(),entry_index=integer(),feature_index=integer(),feature_type=character(),description=character(),source_feature_id=character(),start=integer(),end=integer(),start_modifier=character(),end_modifier=character(),sequence_context=character(),coordinate_status=character(),evidence_codes=character(),evidence_json=character(),location_json=character())
  eb<-data.frame(record_id=integer(),accession=character(),bond_left=integer(),bond_right=integer(),feature_key=character(),feature_type=character(),description=character(),boundary=character(),internal_bond=logical(),evidence_codes=character(),evidence_json=character(),uniprot_accession=character(),source_feature_id=character())
  records<-pm_up_rows(record_rows,er);features<-pm_up_rows(feature_rows,ef);boundaries<-pm_up_rows(boundary_rows,eb)
  index<-split(seq_len(nrow(boundaries)),paste(boundaries$record_id,boundaries$bond_left,sep=":"))
  site_rows<-terminal_rows<-list();si<-ti<-0L
  sites<-result$cleavage$sites
  if(!is.null(sites)) for(i in seq_len(nrow(sites))) {
    idx<-index[[paste(sites$record_id[i],sites$bond_left[i],sep=":")]]
    if(length(idx)) {
      rows<-boundaries[idx,,drop=FALSE];rows<-rows[rows$internal_bond,,drop=FALSE]
      if(nrow(rows)){rows$site_id<-rep(sites$site_id[i],nrow(rows));si<-si+1L;site_rows[[si]]<-rows}
    }
  }
  for(i in seq_len(nrow(m))) for(role in c("N","C")) {
    k<-if(role=="N") m$start[i]-1L else m$end[i]
    idx<-index[[paste(m$record_id[i],k,sep=":")]]
    if(length(idx)) {
      rows<-boundaries[idx,,drop=FALSE]
      rows$placement_id<-rep(sprintf("placement_%08d",i),nrow(rows));rows$peptide_id<-rep(m$peptide_id[i],nrow(rows))
      rows$role<-rep(role,nrow(rows));rows$start<-rep(m$start[i],nrow(rows));rows$end<-rep(m$end[i],nrow(rows))
      ti<-ti+1L;terminal_rows[[ti]]<-rows
    }
  }
  es<-eb;es$site_id<-character();et<-eb;et$placement_id<-character();et$peptide_id<-character();et$role<-character();et$start<-integer();et$end<-integer()
  sa<-pm_up_rows(site_rows,es);ta<-pm_up_rows(terminal_rows,et)
  result$uniprot<-list(records=records,features=features,site_annotations=sa,terminal_annotations=ta,
    summary=list(mapped_records=nrow(records),sequence_matched_records=sum(records$status=="matched"),
      record_status_counts=as.list(table(records$status)),processing_features=nrow(features),
      site_annotation_rows=nrow(sa),sites_with_annotations=length(unique(sa$site_id)),
      terminal_annotation_rows=nrow(ta),peptide_observations_with_annotations=length(unique(ta$peptide_id)),
      interpretation="Descriptive reference annotations only. No filtering or candidate reclassification. Missing annotations do not establish absence of processing; evidence codes and source references remain explicit.",
      identity_policy="Exact primary accession and identical full sequence required; no suffix stripping, secondary-accession transfer or sequence-coordinate remapping."),
    provenance=loaded$provenance)
  result
}
