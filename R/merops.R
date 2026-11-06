# Optional offline reference annotation. Never changes observation quantities.
pm_read_merops_snapshot <- function(snapshot_path, provenance_path) {
  if (!requireNamespace("digest", quietly=TRUE) || !requireNamespace("jsonlite", quietly=TRUE))
    stop("Install digest and jsonlite for MEROPS snapshot verification.")
  if (!file.exists(snapshot_path) || !file.exists(provenance_path))
    stop("MEROPS snapshot and provenance files must both exist.")
  provenance <- jsonlite::fromJSON(provenance_path, simplifyVector=FALSE)
  hash <- provenance$derived_sha256
  if (!is.character(hash) || length(hash)!=1L || !grepl("^[0-9a-fA-F]{64}$",hash))
    stop("MEROPS provenance must contain derived_sha256 for the exact snapshot file.")
  actual <- digest::digest(file=snapshot_path, algo="sha256")
  if (!identical(tolower(hash),actual)) stop("MEROPS snapshot SHA-256 does not match provenance.")
  if (!identical(provenance$schema_version,"proteomaster-merops-observed-export-1"))
    stop("Unsupported MEROPS snapshot schema_version.")
  # Mark UTF8 strings without converting to the Windows system codepage. Source
  # chemistry/reference text can contain characters absent from that codepage.
  con <- if (grepl("\\.gz$",snapshot_path,ignore.case=TRUE)) gzfile(snapshot_path,"rt") else
    file(snapshot_path,"rt")
  on.exit(close(con))
  reference <- read.csv(con, colClasses="character", na.strings=NULL, check.names=FALSE,
                        stringsAsFactors=FALSE, comment.char="", encoding="UTF-8")
  required <- c("merops_record_id","protease_id","protease_name","substrate_accession","p1_position",
    "substrate_species","p4_p4prime","canonical_octamer","cleavage_type","method_code_raw",
    "additional_evidence_code_raw","reference_raw","reference_ids")
  if (!all(required %in% names(reference)) || anyDuplicated(names(reference)))
    stop("MEROPS snapshot lacks the required compact-schema columns.")
  if (anyNA(reference[required]) || any(!nzchar(reference$merops_record_id)) || anyDuplicated(reference$merops_record_id))
    stop("MEROPS record IDs must be present and unique; required fields cannot contain NA.")
  if (is.null(provenance$records) || length(provenance$records)!=1L ||
      !is.numeric(provenance$records) || provenance$records!=nrow(reference))
    stop("MEROPS record count does not match provenance.")
  valid_octamer <- grepl("^[ACDEFGHIKLMNPQRSTVWY]{8}$",reference$p4_p4prime)
  if (any(!reference$canonical_octamer %in% c("true","false")) ||
      any((reference$canonical_octamer=="true") != valid_octamer))
    stop("MEROPS canonical_octamer flags are inconsistent with sequence context.")
  if (any(nzchar(reference$p1_position) & !grepl("^-?[0-9]+$",reference$p1_position)))
    stop("MEROPS P1 coordinates must be integers or blank.")
  provenance$verified_snapshot_sha256 <- actual
  provenance$verification <- "Snapshot bytes and record count verified against supplied provenance; provenance is not digitally authenticated."
  list(reference=reference,provenance=provenance)
}

pm_annotate_merops <- function(result, snapshot_path, provenance_path) {
  sites <- result$cleavage$sites
  required <- c("site_id","accession","bond_left","bond_right","context_display","context_complete")
  if (!is.data.frame(sites) || !all(required %in% names(sites)))
    stop("Build cleavage sites before adding MEROPS annotations.")
  if (anyNA(sites[required]) || anyDuplicated(sites$site_id) || !is.logical(sites$context_complete) ||
      any(sites$bond_left<1 | sites$bond_left!=floor(sites$bond_left) | sites$bond_right!=sites$bond_left+1))
    stop("Cleavage sites require unique IDs, valid internal bonds and explicit context completeness.")
  loaded <- pm_read_merops_snapshot(snapshot_path,provenance_path)
  ref <- loaded$reference
  # Exact accession identity is deliberate: isoforms are not collapsed.
  p1 <- suppressWarnings(as.numeric(ref$p1_position))
  coordinate_ok <- nzchar(ref$substrate_accession) & !is.na(p1) & p1>0
  coord_index <- split(which(coordinate_ok), paste(ref$substrate_accession[coordinate_ok],p1[coordinate_ok],sep="\t"))
  motif_index <- split(which(ref$canonical_octamer=="true"),ref$p4_p4prime[ref$canonical_octamer=="true"])
  output <- list(); oi <- 0L
  for (i in seq_len(nrow(sites))) {
    site <- sites[i,,drop=FALSE]
    query <- gsub("|","",site$context_display,fixed=TRUE)
    complete <- isTRUE(site$context_complete) && grepl("^[ACDEFGHIKLMNPQRSTVWY]{8}$",query) &&
      grepl("^[ACDEFGHIKLMNPQRSTVWY]{4}\\|[ACDEFGHIKLMNPQRSTVWY]{4}$",site$context_display)
    coord <- coord_index[[paste(site$accession,site$bond_left,sep="\t")]]
    motif <- if (complete) motif_index[[query]] else integer()
    # A record matching both criteria appears once as a reference-site candidate.
    motif <- setdiff(motif,coord)
    idx <- c(coord,motif)
    if (!length(idx)) next
    relation <- c(rep("reference_site_candidate",length(coord)),rep("motif_compatibility",length(motif)))
    status <- rep("unknown",length(idx))
    comparable <- complete & ref$canonical_octamer[idx]=="true"
    status[comparable] <- ifelse(ref$p4_p4prime[idx[comparable]]==query,"agreement","mismatch")
    oi <- oi+1L
    output[[oi]] <- cbind(data.frame(site_id=site$site_id,site_accession=site$accession,
      site_bond_left=site$bond_left,site_bond_right=site$bond_right,site_context=site$context_display,
      relation_type=relation,context_status=status,stringsAsFactors=FALSE),ref[idx,,drop=FALSE])
  }
  empty <- data.frame(site_id=character(),site_accession=character(),site_bond_left=integer(),
    site_bond_right=integer(),site_context=character(),relation_type=character(),context_status=character())
  matches <- if (length(output)) do.call(rbind,output) else cbind(empty,ref[FALSE,,drop=FALSE])
  rownames(matches) <- NULL
  result$merops <- list(matches=matches,summary=list(sites_examined=nrow(sites),reference_records=nrow(ref),
    annotation_rows=nrow(matches),sites_with_annotations=length(unique(matches$site_id)),
    reference_site_candidates=sum(matches$relation_type=="reference_site_candidate"),
    reference_site_context_agreements=sum(matches$relation_type=="reference_site_candidate" & matches$context_status=="agreement"),
    reference_site_context_mismatches=sum(matches$relation_type=="reference_site_candidate" & matches$context_status=="mismatch"),
    reference_site_context_unknown=sum(matches$relation_type=="reference_site_candidate" & matches$context_status=="unknown"),
    motif_compatibility_rows=sum(matches$relation_type=="motif_compatibility"),
    interpretation="Reference annotation and motif compatibility do not identify an active protease or prove cleavage in this sample. Source cleavage type and evidence are retained; no annotation is not evidence of no cleavage."),
    provenance=loaded$provenance)
  result
}
