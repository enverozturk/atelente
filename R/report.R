pm_write_csv <- function(x, path) {
  # write.table otherwise escapes non-native characters before fileEncoding applies.
  if (!isTRUE(l10n_info()[["UTF-8"]])) {
    original <- Sys.getlocale("LC_CTYPE")
    on.exit(Sys.setlocale("LC_CTYPE", original), add = TRUE)
    candidates <- if (.Platform$OS.type == "windows") ".UTF-8" else c("C.UTF-8", "en_US.UTF-8", "UTF-8")
    for (locale in candidates) {
      if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", locale)))) break
    }
    if (!isTRUE(l10n_info()[["UTF-8"]])) stop("A UTF-8 locale is required for CSV export.")
  }
  # UTF-8 BOM for Excel; values remain unaggregated and missing cells remain blank.
  con <- file(path, open = "wb")
  writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  close(con)
  suppressWarnings(write.table(x, path, append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE,
                              na = "", qmethod = "double", fileEncoding = "UTF-8"))
}
pm_report_proteins <- function(proteins, source_proteins = NULL) {
  proteins$fasta_gene <- vapply(proteins$header, function(h) {
    m <- regmatches(h, regexec("(?:^| )GN=([^ ]+)", h, perl = TRUE))[[1]]
    if (length(m)) m[2] else ""
  }, character(1))
  proteins$pd_gene <- rep("", nrow(proteins))
  if (!is.null(source_proteins) && all(c("Accession", "Gene Symbol") %in% names(source_proteins))) {
    valid <- !is.na(source_proteins$Accession) & !is.na(source_proteins[["Gene Symbol"]])
    lookup <- split(as.character(source_proteins[["Gene Symbol"]][valid]), as.character(source_proteins$Accession[valid]))
    proteins$pd_gene <- vapply(proteins$accession, function(a) paste(unique(lookup[[a]]), collapse = "; "), character(1))
  }
  proteins$gene_symbol <- vapply(seq_len(nrow(proteins)), function(i) {
    genes <- trimws(unlist(strsplit(c(proteins$fasta_gene[i], proteins$pd_gene[i]), ";", fixed = TRUE)))
    paste(unique(genes[nzchar(genes)]), collapse = "; ")
  }, character(1))
  proteins
}
pm_report_peptides <- function(peptides, source_peptides) {
  fields <- c("peptide_id", "source_row", "sequence", "modifications", "modifications_available",
              "mapping_status", "match_count", "exclusion_reason", "contaminant_screening")
  out <- peptides[, fields, drop = FALSE]
  out$pd_master_accessions <- if ("Master Protein Accessions" %in% names(source_peptides))
    as.character(source_peptides[["Master Protein Accessions"]]) else rep(NA_character_, nrow(out))
  if ("Leading razor protein" %in% names(source_peptides))
    out$leading_razor_accessions <- as.character(source_peptides[["Leading razor protein"]])
  out
}
pm_report_intensities <- function(measurements) {
  if (is.null(measurements) || !nrow(measurements)) return(data.frame())
  fields <- c("peptide_id","entity_id","comparison","analysis_condition")
  take <- measurements$usable_pair %in% TRUE & is.finite(measurements$peptide_log2_value) &
    is.finite(measurements$protein_log2_value) & complete.cases(measurements[fields])
  if (!any(take)) return(data.frame())
  means <- aggregate(measurements[take,c("peptide_log2_value","protein_log2_value"),drop=FALSE],
                     measurements[take,fields,drop=FALSE], mean)
  linear_mean <- function(x) {
    y <- 2^x
    y[!is.finite(y) | y == 0] <- NA_real_
    y
  }
  means$peptide_intensity <- linear_mean(means$peptide_log2_value)
  means$protein_intensity <- linear_mean(means$protein_log2_value)
  means[c(fields,"peptide_intensity","protein_intensity")]
}
pm_report_selection <- function(pair_view) {
  if (is.null(pair_view)) return(NULL)
  comparisons <- pair_view$comparisons
  if (is.null(comparisons) || !nrow(comparisons)) return(pair_view)
  configured <- pair_view$selection_k
  if (!is.null(configured) && (!is.numeric(configured) || length(configured) != 1L ||
      is.na(configured) || !is.finite(configured) || configured < 1 || configured != floor(configured)))
    stop("representative_selection_k must be a positive integer.")
  stored <- pair_view[["selection", exact=TRUE]]$k_by_comparison
  samples <- pair_view$samples
  ks <- setNames(lapply(seq_len(nrow(comparisons)), function(i) {
    name <- comparisons$Name[i]
    value <- if (!is.null(configured)) configured else stored[[name]]
    if (is.null(value)) {
      counts <- vapply(c(comparisons$Control[i], comparisons$Treatment[i]), function(condition)
        sum(samples$Use == "Yes" & samples$Condition == condition, na.rm=TRUE), integer(1))
      value <- max(1L, min(3L, counts))
    }
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value))
      stop("Invalid stored representative selection k.")
    as.integer(value)
  }), comparisons$Name)
  pair_view$selection <- list(policy="ends_coverage_overlap_intensity_v3", k_by_comparison=ks,
    overlap="Distinct other accepted bonds crossed on the same FASTA record; focal bond excluded",
    intensity="Pooled joint log2 geometric intensity over included common positive samples",
    fallback="No eligible pair: historical coverage and stable-ID ranking",
    inference="Original pair tests and BH families unchanged")
  pair_view$selection_k <- NULL
  pair_view
}
pm_report_pack_tables <- function(x) {
  if (is.data.frame(x)) {
    columns <- lapply(x, function(v) {
      if (is.factor(v)) v <- as.character(v)
      type <- if (is.integer(v)) "integer" else if (is.numeric(v)) "numeric" else if (is.logical(v)) "logical" else if (is.character(v)) "character" else stop("Unsupported report column type.")
      if (is.character(v) && length(v) > 32L) {
        dictionary <- unique(v[!is.na(v)])
        if (length(dictionary) < length(v)/2)
          return(list(type=type, dictionary=I(dictionary), values=I(match(v,dictionary)-1L)))
      }
      list(type=type, values=I(v))
    })
    return(list(`__atelente_table__`=1L, nrow=nrow(x), columns=columns))
  }
  if (is.list(x)) return(lapply(x, pm_report_pack_tables))
  x
}
pm_report_unpack_tables <- function(x) {
  if (is.list(x) && identical(x[["__atelente_table__"]],1L)) {
    n <- x$nrow
    if (!is.numeric(n)||length(n)!=1L||!is.finite(n)||n<0||n!=floor(n)||!is.list(x$columns)) stop("Invalid report table.")
    columns <- lapply(x$columns, function(column) {
      if(!is.list(column)||!is.character(column$type)||length(column$type)!=1L||!column$type %in% c("character","logical","integer","numeric"))stop("Invalid report column type.")
      dictionary <- column$dictionary
      code <- !is.null(dictionary)
      if(code && column$type!="character")stop("Invalid report dictionary type.")
      values <- column$values
      template <- if(code)NA_real_ else switch(column$type,character=NA_character_,logical=NA,integer=NA_integer_,numeric=NA_real_)
      if(is.list(values)) values <- vapply(values,function(v)if(is.null(v))template else v,template)
      if(code){
        if(is.list(dictionary))dictionary <- vapply(dictionary,function(v){if(!is.character(v)||length(v)!=1L||is.na(v))stop("Invalid report dictionary.");v},character(1))
        if(is.logical(values)&&all(is.na(values)))values<-as.numeric(values)
        if(!is.character(dictionary) || !is.numeric(values))stop("Invalid report dictionary.")
        if(any(!is.na(values)&(!is.finite(values)|values<0|values>=length(dictionary)|values!=floor(values)))) stop("Invalid report dictionary index.")
        values <- dictionary[values+1L]
      }
      missing_only <- is.logical(values)&&all(is.na(values))
      valid_type <- switch(column$type,
        character=is.character(values)||missing_only,
        logical=is.logical(values),
        numeric=is.numeric(values)||missing_only,
        integer=(is.numeric(values)||missing_only)&&all(is.na(values)|(is.finite(values)&values==floor(values)&abs(values)<=.Machine$integer.max)))
      if(!valid_type)stop("Report column values do not match their declared type.")
      values <- switch(column$type,character=as.character(values),logical=as.logical(values),integer=as.integer(values),numeric=as.numeric(values))
      if(length(values)!=n)stop("Report column length mismatch.")
      values
    })
    return(structure(columns,class="data.frame",row.names=if(n)seq_len(n) else integer()))
  }
  if (is.list(x) && !is.data.frame(x)) return(lapply(x,pm_report_unpack_tables))
  x
}
pm_report_validate_schema <- function(value) {
  schema <- value$payload_schema
  if (!is.null(schema) && !identical(schema,"atelente-report-data-2")) stop("Unsupported report data schema. Use a compatible Atelente version.")
  if (!is.null(schema)) value <- pm_report_unpack_tables(value)
  value
}
pm_report_encode_json <- function(json, threshold = 1024^2) {
  # Transport only: preserve the complete UTF-8 JSON, including all original rows.
  bytes <- charToRaw(enc2utf8(json))
  if (length(bytes) < threshold) return(json)
  # Despite its R name, memCompress(type="gzip") writes a zlib/deflate stream.
  compressed <- memCompress(bytes, type="gzip")
  envelope <- as.character(jsonlite::toJSON(list(schema="atelente-report-envelope-1",
    encoding="deflate-base64", uncompressed_bytes=length(bytes),
    data=gsub("[[:space:]]", "", jsonlite::base64_enc(compressed))), auto_unbox=TRUE))
  if (nchar(envelope,type="bytes") < length(bytes)) envelope else json
}
pm_report_decode_json <- function(text, ..., max_uncompressed_bytes = Inf) {
  if (!is.numeric(max_uncompressed_bytes)||length(max_uncompressed_bytes)!=1L||
      is.na(max_uncompressed_bytes)||max_uncompressed_bytes<=0)
    stop("Invalid report decode size limit.")
  size_error <- function() stop("Report payload exceeds the refresh decode limit. The source report is unchanged. Re-run from the original inputs or use a separately validated chunked migration; a larger limit is not a memory guarantee.")
  if (nchar(text,type="bytes")>max_uncompressed_bytes) size_error()
  args <- list(...)
  # Structural arrays remain arrays until typed-table markers have been read.
  # Scalar vectors still follow simplifyVector, including explicit FALSE callers.
  # This avoids auto-coercing a list of identically shaped tables into one table.
  structural <- function(json) {
    options <- args; options$simplifyDataFrame <- FALSE; options$simplifyMatrix <- FALSE
    do.call(jsonlite::fromJSON,c(list(txt=json),options))
  }
  value <- structural(text)
  if (identical(value$schema, "atelente-report-envelope-1")) {
    if (!identical(value$encoding,"deflate-base64") || !is.character(value$data) || length(value$data)!=1L ||
        !is.numeric(value$uncompressed_bytes) || length(value$uncompressed_bytes)!=1L ||
        !is.finite(value$uncompressed_bytes) || value$uncompressed_bytes<0)
      stop("Invalid compressed report payload.")
    # The legacy PXD payload is 1.536 GB before its row objects are constructed.
    # Reject oversized declared payloads before decompression/object allocation.
    if(value$uncompressed_bytes>max_uncompressed_bytes) size_error()
    bytes <- memDecompress(jsonlite::base64_dec(value$data), type="gzip")
    if (length(bytes) != value$uncompressed_bytes) stop("Report payload length check failed.")
    text <- rawToChar(bytes)
    Encoding(text) <- "UTF-8"
    value <- structural(text)
  }
  if(is.null(value$payload_schema))return(do.call(jsonlite::fromJSON,c(list(txt=text),args)))
  pm_report_validate_schema(value)
}
pm_report_reference_html <- function(access) {
  if (is.null(access)) return("")
  escape <- function(x) {
    if(is.null(x)||!length(x)||is.na(x[1])||!nzchar(as.character(x[1]))) return("Not declared by source")
    x <- enc2utf8(as.character(x[1]))
    for(pair in list(c("&","&amp;"),c("<","&lt;"),c(">","&gt;"),c('"',"&quot;"),c("'","&#39;"))) x <- gsub(pair[1],pair[2],x,fixed=TRUE)
    x
  }
  if(identical(access$enabled,FALSE)) return('<p class="table-note" role="status">Annotations disabled by user (annotations = FALSE).</p>')
  unsupported <- as.character(unlist(access$unsupported_uniprot_ids,use.names=FALSE))
  note <- if(length(unsupported)) paste0('<p class="table-note">UniProt lookup unavailable for ',length(unsupported),
    ' non-UniProt identifier(s): ',escape(paste(head(unsupported,5L),collapse=", ")),
    if(length(unsupported)>5L) " ..." else "",'. Mappings are retained; full identifiers are recorded in provenance.</p>') else ""
  rows <- vapply(access$providers,function(p) paste0("<tr>",paste0("<td style=\"overflow-wrap:anywhere;white-space:normal\">",
    vapply(list(toupper(p$provider),p$access,p$release,p$source_date,p$retrieved_utc,p$snapshot_sha256),escape,character(1)),"</td>",collapse=""),"</tr>"),character(1))
  paste0('<section class="card" id="referenceAccess"><details><summary>Reference sources</summary><div class="scroll"><table><thead><tr>',
    paste0("<th>",c("Source","Access","Release","Source date","Retrieved (UTC)","Snapshot SHA-256"),"</th>",collapse=""),
    "</tr></thead><tbody>",paste(rows,collapse=""),"</tbody></table></div>",note,"</details></section>")
}
# Canonical report summaries. Measurements are summarized once per directed form
# pair and comparison; mapping-specific preference and selection remain separate.
pm_report_pair_summaries <- function(payload) {
  pv <- payload$pair_view
  required <- c("site_id","record_id","accession","bond_left","bond_right","role",
    "neo_peptide_id","spanning_peptide_id","comparison","control","treatment",
    "neo_log2_fc","spanning_log2_fc","log2_fc","control_common","control_total",
    "treatment_common","treatment_total","status","coverage_pass","selection_k",
    "eligible","preferred","overlap","intensity","minco","totalco",
    "assignment_eligible","Representative","selection_rank","selection_fallback",
    "arithmetic_valid","neo_all_log2_fc","spanning_all_log2_fc","fold_change")
  condition_fields <- c("neo","spanning","neoMissing","spanningMissing","neoZero","spanningZero",
    "meanNeoLog2","meanSpanLog2","meanLog2Ratio","meanAllNeoLog2","meanAllSpanLog2")
  required <- c(required, paste0("control_",condition_fields),paste0("treatment_",condition_fields))
  key <- function(x, cols) {
    if (!is.data.frame(x) || !nrow(x)) return(character())
    do.call(paste,c(x[cols],sep="\034"))
  }
  paircols <- c("site_id","role","neo_peptide_id","spanning_peptide_id")
  cached <- pv$prepared_summaries
  if (!is.null(cached)) {
    if (!is.data.frame(cached)) stop("Prepared pair summaries must be a table.")
    if (!nrow(cached)) return(cached)
    if (!all(required %in% names(cached))) stop("Incomplete prepared pair summary schema.")
    if (anyDuplicated(key(cached,c(paircols,"comparison")))) stop("Duplicate prepared pair summaries.")
    if (anyNA(cached$assignment_eligible) || anyNA(cached$Representative) ||
        any(cached$Representative & !cached$assignment_eligible)) stop("Invalid prepared representative assignment.")
    available <- pv$comparisons$Name
    if (any(!is.na(cached$comparison) & !cached$comparison %in% available)) stop("Unknown comparison in prepared pair summaries.")
    for (name in unique(cached$comparison[!is.na(cached$comparison)])) {
      cp <- pv$comparisons[match(name,pv$comparisons$Name),,drop=FALSE];ii <- which(cached$comparison==name)
      if (any(cached$control[ii]!=cp$Control | cached$treatment[ii]!=cp$Treatment)) stop("Prepared comparison labels differ from the sample design.")
      ks <- pv$selection$k_by_comparison[[name]]
      if (!is.null(ks) && any(cached$selection_k[ii]!=ks)) stop("Prepared selection coverage differs from the frozen analysis setting.")
    }
    for (ii in split(which(cached$assignment_eligible),key(cached[cached$assignment_eligible,,drop=FALSE],c("site_id","role","comparison")))) {
      ranks <- cached$selection_rank[ii]
      if (anyNA(ranks) || !identical(sort(as.integer(ranks)),seq_along(ii)) ||
          sum(cached$Representative[ii])!=1L || !cached$Representative[ii[which(ranks==1L)]]) stop("Invalid prepared representative ranking.")
    }
    return(cached)
  }
  ratios <- pv$ratios
  if (is.null(ratios) || !nrow(ratios)) return(data.frame())
  no_comparison <- is.null(pv$comparisons) || !nrow(pv$comparisons)
  comparisons <- if(no_comparison)data.frame(Name=NA_character_,Control="",Treatment="",stringsAsFactors=FALSE) else pv$comparisons
  make_index <- function(x,cols) {
    z <- split(seq_len(nrow(x)),key(x,cols))
    if (length(z)) list2env(z,hash=TRUE,parent=emptyenv()) else new.env(hash=TRUE,parent=emptyenv())
  }
  get_index <- function(env,k) get0(k,envir=env,inherits=FALSE,ifnotfound=integer())
  formcols <- c("neo_peptide_id","spanning_peptide_id")
  pkeys <- key(ratios,paircols)
  first <- which(!duplicated(pkeys))
  map <- ratios[first,c("site_id","record_id","accession","bond_left","bond_right","role",formcols),drop=FALSE]
  mapkeys <- pkeys[first]
  map_order <- order(mapkeys)
  map <- map[map_order,,drop=FALSE]; mapkeys <- mapkeys[map_order]
  included <- pv$samples$Use=="Yes" & (no_comparison | pv$samples$Condition %in% c(comparisons$Control,comparisons$Treatment))
  included_columns <- unique(pv$samples[["Peptide column"]][included])
  included_counts <- tabulate(match(pkeys[ratios$source_column %in% included_columns],mapkeys),nbins=length(mapkeys))
  if(any(included_counts!=length(included_columns)))stop("Missing or duplicate included sample rows for a mapped pair.")
  # A mapping cannot change the underlying form/sample measurement. Verify this
  # before deduplicating, rather than hiding an inconsistent mapping in a cache.
  samplekeys <- key(ratios,c(formcols,"source_column"))
  canonical <- match(samplekeys,samplekeys)
  observed_fields <- c("neo_status","spanning_status","neo_value","spanning_value","log2_ratio")
  for (field in observed_fields) {
    x <- ratios[[field]]; y <- x[canonical]
    same <- (is.na(x)&is.na(y)) | (!is.na(x)&!is.na(y)&x==y)
    if (!all(same)) stop("Inconsistent quantities for a mapped form/sample pair.")
  }
  if (anyDuplicated(paste(pkeys,ratios$source_column,sep="\034"))) stop("Duplicate pair/sample rows in workbook input.")
  unique_rows <- which(!duplicated(samplekeys))
  rr <- ratios[unique_rows,c(formcols,"source_column",observed_fields),drop=FALSE]
  rm(samplekeys,canonical,pkeys)
  formkeys <- key(rr,formcols)
  ri <- split(seq_len(nrow(rr)),formkeys)
  map_form <- match(key(map,formcols),names(ri))
  links <- pv$links
  li <- match(mapkeys,key(links,paircols))
  map$assignment_eligible <- if (is.null(links)) rep(TRUE,nrow(map)) else
    !is.na(li) & links$assignment_eligible[li] %in% TRUE
  e <- payload$cleavage$evidence; terms <- payload$cleavage$termini; sites <- payload$cleavage$sites
  if (is.null(e)) e <- data.frame()
  if (is.null(terms)) terms <- data.frame()
  e <- e[e$role=="spanning" | e$is_candidate %in% TRUE,,drop=FALSE]
  ei <- make_index(e,c("site_id","record_id","role","peptide_id"))
  ti <- make_index(terms,c("record_id","peptide_id","start","end","role"))
  profile <- payload$cleavage$profile
  known <- identical(profile$enzyme,"trypsin") && is.logical(profile$lysine_blocked) && length(profile$lysine_blocked)==1
  # Determine compatibility once per evidence placement and terminal role.
  compatible <- function(role) vapply(seq_len(nrow(e)),function(j) {
    if (!known) return(FALSE)
    kk <- paste(e$record_id[j],e$peptide_id[j],e$start[j],e$end[j],role,sep="\034")
    ii <- get_index(ti,kk)
    length(ii)>0 && all(!terms$processing_status[ii] %in% c("reference_end","possible_initiator_met_removal") &
      !(role=="N" & terms$pre_digest_n_label[ii] %in% TRUE) & terms$digestion_compatible[ii] %in% TRUE)
  },logical(1))
  cn <- compatible("N"); cc <- compatible("C")
  bonds <- split(sites$bond_left,sites$record_id)
  bonds <- if(length(bonds))list2env(bonds,hash=TRUE,parent=emptyenv()) else new.env(parent=emptyenv())
  preferred <- logical(nrow(map)); overlap <- rep(Inf,nrow(map))
  for (j in seq_len(nrow(map))) {
    ek <- function(role,id) paste(map$site_id[j],map$record_id[j],role,id,sep="\034")
    ni <- get_index(ei,ek(map$role[j],map$neo_peptide_id[j]))
    si <- get_index(ei,ek("spanning",map$spanning_peptide_id[j]))
    preferred[j] <- length(ni)>0 && length(si)>0 && all(if(map$role[j]=="N")cc[ni] else cn[ni]) && all(cn[si]) && all(cc[si])
    focal <- map$bond_left[j]
    if (length(si) && all(e$start[si]<=focal & e$end[si]>focal)) {
      bb <- get_index(bonds,as.character(map$record_id[j]))
      overlap[j] <- length(unique(unlist(lapply(si,function(i)bb[bb>=e$start[i]&bb<e$end[i]&bb!=focal]))))
    }
  }
  map$preferred <- preferred; map$overlap <- overlap
  outputs <- vector("list",nrow(comparisons))
  for (ci in seq_len(nrow(comparisons))) {
    cp <- comparisons[ci,,drop=FALSE]
    samples <- pv$samples[pv$samples$Use=="Yes" & (no_comparison | pv$samples$Condition %in% c(cp$Control,cp$Treatment)),,drop=FALSE]
    if(anyDuplicated(samples[["Peptide column"]]))stop("Duplicate included sample columns.")
    ks <- if(no_comparison)0L else pv$selection$k_by_comparison[[cp$Name]]
    display_k <- max(1,min(3,sum(samples$Condition==cp$Control),sum(samples$Condition==cp$Treatment)))
    if (is.null(ks)) ks <- display_k
    summaries <- lapply(ri,function(ii) {
      x <- rr[ii,,drop=FALSE]
      matched <- match(samples[["Peptide column"]],x$source_column)
      if(anyNA(matched))stop("Missing included sample row in pair summaries.")
      x <- x[matched,,drop=FALSE]
      observed <- function(condition) {
        z <- x[if(no_comparison)FALSE else samples$Condition==condition,,drop=FALSE]
        n <- z$neo_status %in% "observed" & is.finite(z$neo_value) & z$neo_value>0
        s <- z$spanning_status %in% "observed" & is.finite(z$spanning_value) & z$spanning_value>0
        co <- n&s&is.finite(z$log2_ratio)
        mn <- function(v) if(length(v))mean(v) else NA_real_
        c(total=nrow(z),common=sum(co),neo=sum(n),spanning=sum(s),
          neoMissing=sum(z$neo_status %in% "missing"),spanningMissing=sum(z$spanning_status %in% "missing"),
          neoZero=sum(z$neo_status %in% "zero"),spanningZero=sum(z$spanning_status %in% "zero"),
          meanNeoLog2=mn(log2(z$neo_value[co])),meanSpanLog2=mn(log2(z$spanning_value[co])),meanLog2Ratio=mn(z$log2_ratio[co]),
          meanAllNeoLog2=mn(log2(z$neo_value[n])),meanAllSpanLog2=mn(log2(z$spanning_value[s])))
      }
      a <- observed(cp$Control); b <- observed(cp$Treatment)
      positive <- x$neo_status %in% "observed" & x$spanning_status %in% "observed" & is.finite(x$neo_value) & is.finite(x$spanning_value) & x$neo_value>0 & x$spanning_value>0 & is.finite(x$log2_ratio)
      both <- x$neo_status %in% "observed" & x$spanning_status %in% "observed" & is.finite(x$neo_value) & is.finite(x$spanning_value) & x$neo_value>0 & x$spanning_value>0
      if(any(both & !is.finite(x$log2_ratio)) || any(!both & is.finite(x$log2_ratio)))stop("Pair ratio is inconsistent with its measured components.")
      expected <- log2(x$neo_value[positive])-log2(x$spanning_value[positive])
      if(any(abs(expected-x$log2_ratio[positive])>1e-10*pmax(1,abs(expected),abs(x$log2_ratio[positive]))))stop("Inconsistent sample ratio in workbook input.")
      c(setNames(a,paste0("control_",names(a))),setNames(b,paste0("treatment_",names(b))),
        neo_log2_fc=unname(b["meanNeoLog2"]-a["meanNeoLog2"]),spanning_log2_fc=unname(b["meanSpanLog2"]-a["meanSpanLog2"]),
        log2_fc=unname(b["meanLog2Ratio"]-a["meanLog2Ratio"]),neo_all_log2_fc=unname(b["meanAllNeoLog2"]-a["meanAllNeoLog2"]),
        spanning_all_log2_fc=unname(b["meanAllSpanLog2"]-a["meanAllSpanLog2"]),
        intensity=if(any(positive))mean((log2(x$neo_value[positive])+log2(x$spanning_value[positive]))/2) else -Inf)
    })
    numeric <- as.data.frame(do.call(rbind,summaries),stringsAsFactors=FALSE)
    rownames(numeric) <- NULL
    x <- cbind(map,numeric[map_form,,drop=FALSE]); rownames(x) <- NULL
    x$comparison <- cp$Name; x$control <- if(no_comparison)NA_character_ else cp$Control; x$treatment <- if(no_comparison)NA_character_ else cp$Treatment
    x$selection_k <- ks; x$arithmetic_valid <- TRUE
    # Missing and zero remain distinct. Display k is independent of the frozen
    # analysis k and cannot replace the representative selected below.
    fc <- is.finite(x$log2_fc)
    both_t <- !fc & x$treatment_common>0 & x$control_total>0 & x$control_neoMissing==x$control_total & x$control_spanningMissing==x$control_total
    both_c <- !fc & !both_t & x$control_common>0 & x$treatment_total>0 & x$treatment_neoMissing==x$treatment_total & x$treatment_spanningMissing==x$treatment_total
    neo_t <- !fc & !both_t & !both_c & x$treatment_common>0 & x$control_total>0 & x$control_neoMissing==x$control_total & x$control_spanning>0
    neo_c <- !fc & !both_t & !both_c & !neo_t & x$control_common>0 & x$treatment_total>0 & x$treatment_neoMissing==x$treatment_total & x$treatment_spanning>0
    span_only <- !fc & !both_t & !both_c & !neo_t & !neo_c & x$control_neo+x$treatment_neo==0 & x$control_spanning+x$treatment_spanning>0
    neo_only <- !fc & !both_t & !both_c & !neo_t & !neo_c & !span_only & x$control_spanning+x$treatment_spanning==0 & x$control_neo+x$treatment_neo>0
    coverage <- function(k,eligible=FALSE) (fc & x$control_common>=k & x$treatment_common>=k) |
      (both_t & x$treatment_common>=k) | (both_c & x$control_common>=k) |
      (neo_t & x$treatment_common>=k & x$control_spanning>=k) | (neo_c & x$control_common>=k & x$treatment_spanning>=k) |
      (!eligible & span_only & pmax(x$control_spanning,x$treatment_spanning)>=k) |
      (!eligible & neo_only & pmax(x$control_neo,x$treatment_neo)>=k)
    x$status <- "Insufficient data"
    x$status[fc] <- "Quantified in both conditions"
    x$status[both_t] <- paste("Both components quantified only in",cp$Treatment)
    x$status[both_c] <- paste("Both components quantified only in",cp$Control)
    x$status[neo_t] <- paste("Neo quantified only in",cp$Treatment)
    x$status[neo_c] <- paste("Neo quantified only in",cp$Control)
    x$status[span_only] <- "Spanning quantified; neo not quantified"
    x$status[neo_only] <- "Neo quantified; spanning not quantified"
    x$coverage_pass <- coverage(display_k); x$eligible <- coverage(ks,TRUE)
    if(no_comparison){x$status <- "No comparison";x$coverage_pass <- FALSE;x$eligible <- FALSE}
    x$minco <- pmin(x$control_common,x$treatment_common); x$totalco <- x$control_common+x$treatment_common
    x$fold_change <- 2^x$log2_fc; x$fold_change[!is.finite(x$fold_change)|x$fold_change<=0] <- NA_real_
    x$Representative <- FALSE; x$selection_rank <- NA_integer_; x$selection_fallback <- FALSE
    # Same stable group/identity ordering as the original workbook implementation.
    groups <- split(which(x$assignment_eligible),key(x[x$assignment_eligible,,drop=FALSE],c("site_id","role","comparison")))
    ordered <- integer()
    for (ii in groups) {
      fallback <- !any(x$eligible[ii])
      o <- if(fallback)order(-x$minco[ii],-x$totalco[ii],x$neo_peptide_id[ii],x$spanning_peptide_id[ii],method="radix") else
        order(-as.integer(x$eligible[ii]),ifelse(x$eligible[ii],-as.integer(x$preferred[ii]),0),-x$minco[ii],-x$totalco[ii],
          ifelse(x$eligible[ii],x$overlap[ii],0),ifelse(x$eligible[ii],-x$intensity[ii],0),x$neo_peptide_id[ii],x$spanning_peptide_id[ii],method="radix")
      x$Representative[ii[o[1]]] <- TRUE; x$selection_rank[ii[o]] <- seq_along(ii); x$selection_fallback[ii] <- fallback
      ordered <- c(ordered,ii)
    }
    outputs[[ci]] <- x[c(ordered,which(!x$assignment_eligible)),,drop=FALSE]
  }
  out <- do.call(rbind,outputs); rownames(out) <- NULL
  contrasts <- pv$contrasts
  if(is.data.frame(contrasts) && nrow(contrasts) && !no_comparison){
    ckeys <- key(contrasts,c(paircols,"comparison"))
    if(anyDuplicated(ckeys))stop("Duplicate exported pair contrast.")
    j <- match(key(out,c(paircols,"comparison")),ckeys);ii <- which(!is.na(j));jj <- j[ii]
    close <- function(a,b)(is.na(a)&is.na(b)) | (is.finite(a)&is.finite(b)&abs(a-b)<=1e-10*pmax(1,abs(a),abs(b)))
    if("log2_fc" %in% names(contrasts) && !all(close(out$log2_fc[ii],contrasts$log2_fc[jj])))stop("Exported pair effect differs from its sample ratios.")
    if("fold_change" %in% names(contrasts)){
      use <- which(is.finite(contrasts$fold_change[jj]))
      if(length(use) && !all(close(out$log2_fc[ii[use]],log2(contrasts$fold_change[jj[use]]))))stop("Exported pair fold change differs from its sample ratios.")
    }
    suffix <- c(samples="total",neo="neo",spanning="spanning",usable_pairs="common")
    for(side in c("control","treatment"))for(field in names(suffix)){
      old <- paste0(side,"_",field);new <- paste0(side,"_",suffix[[field]])
      if(old %in% names(contrasts) && any(is.na(contrasts[[old]][jj]) | out[[new]][ii]!=contrasts[[old]][jj]))stop("Exported pair coverage differs from its sample design.")
    }
  }
  out
}


pm_write_report <- function(payload, template, path) {
  if (!is.null(payload$pair_view)) payload$pair_view <- pm_report_selection(payload$pair_view)
  if (!is.null(payload$cleavage) && basename(template) == "report_template.html") {
    template <- file.path(dirname(template), "cleavage_report_template.html")
    if (!is.null(payload$quantities)) {
      # A compact attached report already stores this summary. Refresh must
      # preserve it when the original sample rows are no longer embedded.
      intensity_summary <- if (!is.null(payload$quantities$measurements))
        pm_report_intensities(payload$quantities$measurements) else payload$quantities$intensity_summary
      payload$quantities <- payload$quantities[intersect(c("results", "comparisons","paired_bh_counts"), names(payload$quantities))]
      payload$quantities$intensity_summary <- intensity_summary
    }
  }
  if (!is.null(payload$pair_view) && !is.null(payload$cleavage)) {
    payload$pair_view$prepared_summaries <- pm_report_pair_summaries(payload)
    payload$pair_view$ratios <- NULL
    payload$pair_view$contrasts <- NULL
    # Mapping relationships and original test values stay in the report; large
    # modelling/assignment audit columns remain in the separate audit workbook.
    pair_keys <- c("site_id","role","neo_peptide_id","spanning_peptide_id")
    for (table in c("links","statistics")) {
      fields <- if(table=="links")c(pair_keys,"assignment_eligible","exclusion_reason") else
        c(pair_keys,"comparison","family","status","p_value","adjusted_p_value")
      if(is.data.frame(payload$pair_view[[table]]))
        payload$pair_view[[table]] <- payload$pair_view[[table]][intersect(fields,names(payload$pair_view[[table]]))]
    }
    excluded <- payload$pair_view$excluded_statistics
    if (is.data.frame(excluded)) payload$pair_view$excluded_statistics <- excluded[intersect(c("site_id","role","neo_peptide_id","spanning_peptide_id","comparison","reason","family"),names(excluded))]
  }
  payload$payload_schema <- "atelente-report-data-2"
  # JSON readers do not retain AsIs on singleton arrays across refresh.
  if (!is.null(payload$export_workbook)) {
    for (field in c("sheets","internal_sheets"))
      if (!is.null(payload$export_workbook[[field]]))
        payload$export_workbook[[field]] <- I(unname(as.character(payload$export_workbook[[field]])))
  }
  # AsIs significant digits is explicit across jsonlite 1.x and 2.x. digits=NA
  # can use only 15 significant digits in some older builds, losing precision.
  json <- as.character(jsonlite::toJSON(pm_report_pack_tables(payload), dataframe = "rows", auto_unbox = TRUE, digits = I(17), na = "null", null = "null"))
  json <- gsub("<", "\\u003c", json, fixed = TRUE)
  json <- gsub("&", "\\u0026", json, fixed = TRUE)
  # A hard transport guard, not a promise about a browser's available memory.
  # Refuse a single decoded string close to the V8 limit before writing HTML.
  if (nchar(json,type="bytes") > 450*1024^2)
    stop("Report data exceeds the single-file transport limit. Excel/audit outputs are retained; HTML was not written.")
  html <- paste(readLines(template, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  workbook_marker <- "__ATELENTE_WORKBOOK_JS__"
  if (grepl(workbook_marker, html, fixed=TRUE)) {
    workbook_parts <- strsplit(html, workbook_marker, fixed=TRUE)[[1]]
    if (length(workbook_parts)!=2L) stop("Expected exactly one workbook script placeholder.")
    workbook_helper <- file.path(dirname(template), "report_workbook.js")
    if (!file.exists(workbook_helper)) stop("Missing report asset: report_workbook.js")
    html <- paste0(workbook_parts[1],paste(readLines(workbook_helper,warn=FALSE,encoding="UTF-8"),collapse="\n"),workbook_parts[2])
  }
  payload_marker <- "__ATELENTE_REPORT_PAYLOAD_JS__"
  if (grepl(payload_marker, html, fixed=TRUE)) {
    payload_parts <- strsplit(html, payload_marker, fixed=TRUE)[[1]]
    if (length(payload_parts)!=2L) stop("Expected exactly one payload-codec script placeholder.")
    codec <- file.path(dirname(template), "report_payload.js")
    if (!file.exists(codec)) stop("Missing report asset: report_payload.js")
    html <- paste0(payload_parts[1], paste(readLines(codec,warn=FALSE,encoding="UTF-8"),collapse="\n"), payload_parts[2])
    json <- pm_report_encode_json(json)
  }
  # Keep the report self-contained in both source and installed-package workflows.
  js_marker <- "__ATELENTE_PAIR_SUMMARY_JS__"
  if (grepl(js_marker, html, fixed = TRUE)) {
    js_parts <- strsplit(html, js_marker, fixed = TRUE)[[1]]
    if (length(js_parts) != 2L) stop("Expected exactly one pair-summary script placeholder.")
    helper <- file.path(dirname(template), "pair_summary.js")
    if (!file.exists(helper)) stop("Missing report asset: pair_summary.js")
    js <- paste(readLines(helper, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    html <- paste0(js_parts[1], js, js_parts[2])
  }
  parts <- strsplit(html, "__PROTEOMASTER_DATA__", fixed = TRUE)[[1]]
  if (length(parts) != 2L) stop("Report template must contain exactly one data placeholder.")
  assembled <- paste0(parts[1], json, parts[2])
  reference_html <- pm_report_reference_html(payload$reference_access)
  if(nzchar(reference_html)) {
    marker <- if(identical(payload$reference_access$enabled,FALSE)) "</header>" else "</main>"
    where <- regexpr(marker,assembled,fixed=TRUE)[1]
    if(where<1L) stop("Report template lacks a reference metadata insertion point.")
    if(marker=="</header>") where <- where+nchar(marker)
    # substring() defaults to last=1e6, which truncates large self-contained reports.
    assembled <- paste0(substr(assembled,1,where-1L),reference_html,substring(assembled,where,nchar(assembled)))
  }
  writeLines(assembled, path, useBytes = TRUE)
}
pm_export <- function(result, output_dir, template = system.file("templates", "report_template.html", package = "atelente", mustWork = TRUE), inputs, metadata = list(), excel = FALSE, sample_abundance_type = NULL, neo_design = NULL, representative_selection_k = NULL, compact_bundle = FALSE) {
  if (!requireNamespace("jsonlite", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE))
    stop("Install jsonlite and digest to export reports.")
  if (!is.logical(excel) || length(excel) != 1L || is.na(excel)) stop("excel must be TRUE or FALSE.")
  if (excel && (!exists("pm_write_excel", mode = "function") || !requireNamespace("openxlsx", quietly = TRUE)))
    stop("Source R/excel.R and install openxlsx for Excel export.")
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("Output directory must be empty; choose a new directory to preserve previous results.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  source <- result$source_peptides
  source <- data.frame(peptide_id = result$peptides$peptide_id, source_row = result$peptides$source_row, source, check.names = FALSE)
  tables <- list(peptides = result$peptides, matches = result$matches, pd_coordinate_checks = result$pd_checks,
                 reported_candidates = result$candidates, measurements = result$measurements,
                 excluded_peptides = result$excluded, source_peptides = source)
  if (!is.null(result$source_input)) tables$source_input <- result$source_input
  if (!is.null(result$column_mapping_audit))
    jsonlite::write_json(result$column_mapping_audit, file.path(output_dir,"column_mapping_audit.json"),
      auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
  if (!is.null(result$cleavage)) {
    tables$cleavage_sites <- result$cleavage$sites
    tables$cleavage_evidence <- result$cleavage$evidence
    tables$neo_spanning_ratios <- pm_neo_spanning_ratios(result$cleavage$evidence, result$measurements, sample_abundance_type)
    if (!is.null(neo_design) && nrow(neo_design$comparisons)) {
      tables$neo_spanning_contrasts <- pm_neo_spanning_contrasts(tables$neo_spanning_ratios, neo_design)
      tables$neo_spanning_samples <- neo_design$samples
    }
    tables$peptide_termini <- result$cleavage$termini
    if (!is.null(result$cleavage$placements)) tables$placement_specificity <- result$cleavage$placements
    if (!is.null(result$cleavage$excluded_placements)) tables$excluded_placements <- result$cleavage$excluded_placements
  }
  if (!is.null(neo_design)) tables$neo_spanning_samples <- neo_design$samples
  if (!is.null(result$pair_statistics)) {
    workflow <- result$pair_statistics
    for (name in c("results","measurement_results","links","excluded","membership"))
      if (is.data.frame(workflow[[name]])) tables[[paste0("pair_",name)]] <- workflow[[name]]
    jsonlite::write_json(list(settings=workflow$settings,warnings=workflow$warnings),
      file.path(output_dir,"pair_statistics_settings.json"),auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
  }
  if (!is.null(result$merops)) tables$merops_matches <- result$merops$matches
  if (!is.null(result$uniprot)) {
    for (name in c("records", "features", "site_annotations", "terminal_annotations"))
      tables[[paste0("uniprot_", name)]] <- result$uniprot[[name]]
  }
  if (!is.null(result$source_proteins)) {
    tables$source_proteins <- result$source_proteins
    p <- result$source_proteins
    ids <- if ("Accession" %in% names(p)) as.character(p$Accession) else rep(NA_character_, nrow(p))
    pf <- if ("Contaminant" %in% names(p)) pm_flag(p$Contaminant) else rep(NA, nrow(p))
    tables$protein_screening <- data.frame(source_row = seq_len(nrow(p))+1L, accession = ids,
      contaminant_flag = pf, excluded = pf %in% TRUE | ids %in% result$settings$contaminant_ids)
  }
  for (name in names(tables)) pm_write_csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")))
  if (excel && !compact_bundle) {
    result$analysis_tables <- tables[intersect(c("neo_spanning_contrasts","pair_results","pair_measurement_results",
      "pair_excluded","placement_specificity","excluded_placements"),names(tables))]
    pm_write_excel(result, file.path(output_dir, "proteomaster_results.xlsx"))
  }
  sample_file <- NULL
  if (excel && exists("pm_sample_template", mode="function") && nrow(result$measurements)) {
    selected_type <- sample_abundance_type
    if (is.null(selected_type) && "abundance" %in% result$measurements$abundance_type) selected_type <- "abundance"
    if (!is.null(selected_type)) {
      pm_write_sample_metadata(pm_sample_template(result, selected_type), file.path(output_dir,"sample_metadata.xlsx"))
      sample_file <- "sample_metadata.xlsx"
    }
  }
  # Metadata is not needed for mapping. Conditions are not inferred from sample names.
  columns <- unique(result$measurements[, c("source_column", "abundance_type", "measurement_label"), drop = FALSE])
  sample_cols <- columns[columns$abundance_type != "grouped", , drop = FALSE]
  if (!nrow(sample_cols)) sample_cols <- columns
  sample_cols$sample_id <- rep("", nrow(sample_cols)); sample_cols$condition <- rep("", nrow(sample_cols))
  pm_write_csv(sample_cols, file.path(output_dir, "sample_metadata_template.csv"))
  status <- as.data.frame(table(result$peptides$mapping_status), stringsAsFactors = FALSE)
  names(status) <- c("status", "rows")
  pm_write_csv(status, file.path(output_dir, "mapping_summary.csv"))
  input_manifest <- lapply(inputs, function(path) list(file = normalizePath(path, winslash = "/", mustWork = TRUE),
    sha256 = digest::digest(file = path, algo = "sha256"), bytes = unname(file.info(path)$size)))
  manifest <- list(tool = "Atelente mapping prototype", version = "0.5.0", report_version = "0.9.0", excel_export = excel, sample_metadata_file = sample_file, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    inputs = input_manifest, settings = result$settings, metadata = metadata,
    peptide_rows = nrow(result$peptides), excluded_rows = nrow(result$excluded), fasta_records = nrow(result$fasta),
    duplicate_fasta_accessions = sum(duplicated(result$fasta$accession)),
    limitations = c("Reported scope does not establish uniqueness across the full proteome.",
      "Exact amino-acid matching only; I and L remain distinct.",
      if (is.null(result$cleavage)) "Cleavage classification was not requested." else "Cleavage candidates are provisional; normal processing is not comprehensively annotated. Absence of a supporting peptide is not evidence against cleavage.",
      if (is.null(neo_design) || !nrow(neo_design$comparisons)) "No neo/spanning contrast was requested." else
        if(!isTRUE(result$pair_statistics$settings$test)) "Neo/spanning contrasts describe changes in observed sample ratios. They carry no significance test or enzyme assignment." else
        "Experimental pair tests use explicit independent units. A test is a measurement pair within its contrast/family, not an independent test for each mapped bond.",
      "Optional absent annotations remain unknown. The supplied FASTA is not assumed to be the original search database."))
  if (!is.null(result$merops)) {
    manifest$merops <- result$merops[c("summary", "provenance")]
    jsonlite::write_json(manifest$merops, file.path(output_dir, "merops_provenance.json"), auto_unbox=TRUE, pretty=TRUE, na="null", null="null")
  }
  if (!is.null(result$uniprot)) {
    manifest$uniprot <- result$uniprot[c("summary", "provenance")]
    jsonlite::write_json(manifest$uniprot, file.path(output_dir, "uniprot_provenance.json"), auto_unbox=TRUE, pretty=TRUE, na="null", null="null")
  }
  if (!is.null(result$reference_access)) manifest$reference_access <- result$reference_access
  if (!is.null(result$pair_statistics)) manifest$pair_statistics <- result$pair_statistics$settings
  jsonlite::write_json(manifest, file.path(output_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE, na = "null", null="null")
  # Include only records used by plotted matches; raw input and quantities remain in CSV.
  used_records <- unique(result$matches$record_id)
  proteins <- result$fasta[match(used_records, result$fasta$record_id), c("record_id", "accession", "header", "length"), drop = FALSE]
  proteins <- pm_report_proteins(proteins, result$source_proteins)
  payload <- list(summary = status, scope = result$settings$scope, fasta_records = nrow(result$fasta),
    peptides = pm_report_peptides(result$peptides, result$source_peptides),
    matches = result$matches, proteins = proteins)
  if (!is.null(result$reference_access)) payload$reference_access <- result$reference_access
  if (!is.null(result$column_mapping_audit)) payload$column_mapping_audit <- result$column_mapping_audit
  if (!is.null(result$cleavage)) payload$cleavage <- result$cleavage[c("sites", "evidence", "profile", "termini", "placements")]
  if (!is.null(result$cleavage) || compact_bundle) {
    report_type <- sample_abundance_type
    if (is.null(report_type)) report_type <- if ("abundance" %in% result$measurements$abundance_type) "abundance" else "normalized"
    report_measurements <- result$measurements[result$measurements$abundance_type %in% intersect(report_type, c("abundance", "normalized")),
      c("peptide_id", "source_column", "value", "value_status"), drop = FALSE]
    sample_quantities <- list(measurements = report_measurements,
      samples = if (!is.null(neo_design)) neo_design$samples else NULL,
      comparisons = if (!is.null(neo_design)) neo_design$comparisons else NULL)
    if (!is.null(result$cleavage)) {
      payload$pair_view <- c(list(ratios = tables$neo_spanning_ratios,
        contrasts = tables$neo_spanning_contrasts), sample_quantities)
    } else {
      # Mapping-only compact exports retain samples without inventing pairs,
      # cleavage evidence, or representative-selection metadata.
      payload$sample_quantities <- sample_quantities
    }
  }
  if (!is.null(payload$pair_view)) payload$pair_view$selection_k <- representative_selection_k
  if (!is.null(result$pair_statistics) && !is.null(payload$pair_view)) {
    payload$pair_view$links <- result$pair_statistics$links
    payload$pair_view$statistics <- result$pair_statistics$results
    payload$pair_view$excluded_statistics <- result$pair_statistics$excluded
    payload$pair_view$inference_settings <- result$pair_statistics$settings
  }
  if (!is.null(result$merops)) payload$merops <- result$merops[c("matches", "summary")]
  if (!is.null(result$uniprot)) payload$uniprot <- result$uniprot[c("records", "site_annotations", "terminal_annotations", "summary")]
  if (exists("pm_site_summary", mode = "function") && !is.null(result$cleavage)) payload$site_summary <- pm_site_summary(result)
  if (!is.null(metadata$dataset_label)) payload$dataset_label <- metadata$dataset_label
  if (excel) payload$excel_file <- "proteomaster_results.xlsx"
  if (!is.null(sample_file)) payload$sample_metadata_file <- sample_file
  downloadable <- c(pair_tests="pair_measurement_results",pair_mapping_results="pair_results",
    pair_exclusions="pair_excluded",sample_ratios="neo_spanning_ratios",pair_contrasts="neo_spanning_contrasts",
    specificity="placement_specificity",excluded_placements="excluded_placements",source_input="source_input")
  downloadable <- downloadable[downloadable %in% names(tables)]
  payload$analysis_files <- as.list(setNames(paste0(unname(downloadable),".csv"),names(downloadable)))
  if (compact_bundle) {
    if (!exists("pm_prepare_export_bundle",mode="function") || !exists("pm_finalize_audit",mode="function"))
      stop("Load export_bundle, report_exports and audit_bundle modules for compact exports.")
    payload <- pm_prepare_export_bundle(payload,output_dir,manifest,excel)
  }
  pm_write_report(payload, template, file.path(output_dir, "report.html"))
  if (compact_bundle) {
    manifest <- pm_finalize_audit(output_dir,manifest)
    public_files <- file.path(output_dir,c("report.html","atelente_results.xlsx","atelente_audit.zip"))
    public_files <- public_files[file.exists(public_files)]
    manifest$output_sha256 <- as.list(setNames(vapply(public_files,function(path)
      digest::digest(file=path,algo="sha256"),character(1)),basename(public_files)))
    jsonlite::write_json(manifest,file.path(output_dir,"run_manifest.json"),auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
    return(invisible(manifest))
  }
  files <- list.files(output_dir,full.names=TRUE)
  files <- files[basename(files)!="run_manifest.json" & !dir.exists(files)]
  manifest$output_sha256 <- as.list(setNames(vapply(files,function(path)
    digest::digest(file=path,algo="sha256"),character(1)),basename(files)))
  jsonlite::write_json(manifest,file.path(output_dir,"run_manifest.json"),auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
  invisible(manifest)
}
