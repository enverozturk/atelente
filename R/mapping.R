# Atelente: peptide placement, not protease inference.
# Coordinates are one-based and inclusive. No I/L equivalence is assumed.

pm_upper <- function(x) chartr("abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", x)

pm_bind <- function(rows, empty) {
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
pm_split <- function(x) {
  if (!length(x) || is.na(x) || !nzchar(trimws(x))) return(character())
  unique(trimws(strsplit(x, ";", fixed = TRUE)[[1]]))
}
pm_value <- function(x, name, i) {
  if (!name %in% names(x)) return(NA_character_)
  as.character(x[[name]][i])
}
pm_flag <- function(x) {
  if (is.logical(x)) return(x)
  x <- trimws(as.character(x))
  out <- rep(NA, length(x))
  out[x %in% c("TRUE", "True", "true", "1", "Yes", "yes", "YES", "DO\u011eRU")] <- TRUE
  out[x %in% c("FALSE", "False", "false", "0", "No", "no", "NO", "YANLI\u015e")] <- FALSE
  out
}
pm_read_table <- function(path, sheet = NULL) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "xlsx") {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the readxl package to read XLSX files.")
    available <- readxl::excel_sheets(path)
    if (is.null(sheet)) {
      if (length(available) != 1L) stop("Choose a worksheet explicitly: ", paste(available, collapse = ", "))
      sheet <- available[1]
    }
    # Import source cells as text, as for CSV/TSV. Guessing from early rows can
    # silently discard sparse annotations or mixed PD fields later in a column.
    # Numeric quantities and flags are interpreted explicitly downstream.
    x <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text",
                                        .name_repair = "minimal"), stringsAsFactors = FALSE)
  } else if (ext %in% c("csv", "tsv", "txt")) {
    # R's UTF-8-BOM text connection transcodes to the current locale. Under
    # LC_CTYPE=C it can stop at the first non-ASCII field and return fewer rows.
    if (!isTRUE(l10n_info()[["UTF-8"]])) {
      original <- Sys.getlocale("LC_CTYPE")
      on.exit(Sys.setlocale("LC_CTYPE", original), add = TRUE)
      locales <- if (.Platform$OS.type == "windows") ".UTF-8" else c("C.UTF-8", "en_US.UTF-8", "UTF-8")
      for (locale in locales) {
        if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", locale)))) break
      }
      if (!isTRUE(l10n_info()[["UTF-8"]])) stop("A UTF-8 locale is required to read delimited tables.")
    }
    sep <- if (ext == "csv") "," else "\t"
    if (ext == "csv") {
      # Count separators outside quoted fields in the first logical record.
      # Only the header selects the delimiter; decimal commas in data do not.
      probe <- readLines(path, n=64L, encoding="UTF-8", warn=FALSE)
      probe <- probe[seq_along(probe)>=match(TRUE,nzchar(trimws(probe)),nomatch=length(probe)+1L)]
      quotes <- vapply(gregexpr('"',probe,fixed=TRUE),function(x)if(x[1]<0L)0L else length(x),integer(1))
      header_end <- which(cumsum(quotes) %% 2L == 0L)[1]
      if(length(probe) && is.na(header_end))
        stop("CSV header has an unclosed quoted field or exceeds 64 lines. Use a single header row or XLSX.",call.=FALSE)
      if(length(probe)) probe <- probe[seq_len(header_end)]
      fields <- function(delimiter) {
        con <- textConnection(probe);on.exit(close(con))
        counts <- count.fields(con, sep=delimiter, quote='"', comment.char="", blank.lines.skip=TRUE)
        counts <- counts[!is.na(counts)]
        if(length(counts)) counts[1] else 0L
      }
      counts <- vapply(c(",",";"),fields,integer(1))
      if (all(counts>1L))
        stop("Ambiguous CSV delimiter: the header contains both unquoted commas and semicolons. Use consistently quoted CSV, TSV or XLSX.",call.=FALSE)
      if (counts[[";"]]>1L) sep <- ";"
    }
    x <- read.table(path, header = TRUE, sep = sep, quote = '"', comment.char = "", check.names = FALSE,
                    stringsAsFactors = FALSE, colClasses = "character", na.strings = "", fileEncoding = "UTF-8-BOM",
                    encoding = "UTF-8")
    sheet <- ""
  } else stop("Supported peptide/protein tables: XLSX, CSV, TSV.")
  if (anyDuplicated(names(x)) || any(!nzchar(names(x)))) stop("Column names must be nonempty and unique.")
  attr(x, "source_file") <- normalizePath(path, winslash = "/", mustWork = TRUE)
  attr(x, "source_sheet") <- sheet
  x
}
pm_read_fasta <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt", encoding = "UTF-8") else file(path, "rt", encoding = "UTF-8")
  on.exit(close(con))
  lines <- trimws(readLines(con, warn = FALSE))
  lines <- lines[nzchar(lines)]
  starts <- which(startsWith(lines, ">"))
  if (!length(starts) || starts[1] != 1L) stop("FASTA must start with a header.")
  ends <- c(starts[-1] - 1L, length(lines))
  if (any(ends == starts)) stop("FASTA contains an empty sequence.")
  headers <- substring(lines[starts], 2L)
  tokens <- sub("[[:space:]].*$", "", headers)
  accessions <- vapply(strsplit(tokens, "|", fixed = TRUE), function(z) {
    if (length(z) >= 3L && z[1] %in% c("sp", "tr")) z[2] else z[1]
  }, character(1))
  sequences <- vapply(seq_along(starts), function(i) pm_upper(paste0(lines[(starts[i]+1L):ends[i]], collapse = "")), character(1))
  if (any(!grepl("^[A-Z*]+$", sequences))) stop("FASTA contains unsupported sequence characters.")
  # Keep duplicate accession records: record_id, not accession, is the unique key.
  data.frame(record_id = seq_along(starts), accession = accessions, header = headers,
             sequence = sequences, length = nchar(sequences), stringsAsFactors = FALSE)
}
pm_parse_sequence <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(list(sequence = NA_character_, left = NA_character_, right = NA_character_, status = "missing_sequence"))
  x <- trimws(x)
  m <- regmatches(x, regexec("^\\[([^]]+)\\]\\.([A-Za-z]+)\\.\\[([^]]+)\\]$", x))[[1]]
  if (length(m)) return(list(sequence = pm_upper(m[3]), left = m[2], right = m[4], status = "parsed"))
  if (grepl("^[A-Za-z]+$", x)) return(list(sequence = pm_upper(x), left = NA_character_, right = NA_character_, status = "parsed"))
  list(sequence = NA_character_, left = NA_character_, right = NA_character_, status = "unsupported_sequence_format")
}
# Column-wise forms of pm_split() and pm_parse_sequence(); element i equals the
# scalar function applied to x[i]. Used to avoid per-row regex/trimws calls.
pm_split_all <- function(x) {
  x <- as.character(x)
  out <- rep(list(character()), length(x))
  keep <- !is.na(x) & nzchar(trimws(x))
  if (any(keep)) {
    parts <- strsplit(x[keep], ";", fixed = TRUE)
    sizes <- lengths(parts)
    flat <- trimws(unlist(parts, use.names = FALSE))
    out[keep] <- lapply(split(flat, rep.int(seq_along(parts), sizes)), unique)[as.character(seq_along(parts))]
    # strsplit("", ...) style empty pieces: an element whose parts are all empty
    # still yields unique(trimws(...)) of those pieces, as in pm_split().
    empty <- which(sizes == 0L)
    if (length(empty)) out[which(keep)[empty]] <- list(character())
  }
  out
}
pm_parse_sequence_all <- function(x) {
  x <- as.character(x); n <- length(x)
  sequence <- left <- right <- rep(NA_character_, n)
  status <- rep("unsupported_sequence_format", n)
  missing <- is.na(x) | !nzchar(trimws(x))
  status[missing] <- "missing_sequence"
  t <- trimws(x)
  m <- regmatches(t, regexec("^\\[([^]]+)\\]\\.([A-Za-z]+)\\.\\[([^]]+)\\]$", t))
  flank <- !missing & lengths(m) > 0L
  if (any(flank)) {
    mm <- do.call(rbind, m[flank])
    sequence[flank] <- pm_upper(mm[, 3]); left[flank] <- mm[, 2]; right[flank] <- mm[, 4]
    status[flank] <- "parsed"
  }
  plain <- !missing & !flank & grepl("^[A-Za-z]+$", t)
  sequence[plain] <- pm_upper(t[plain]); status[plain] <- "parsed"
  list(sequence = sequence, left = left, right = right, status = status)
}
# Parse PD "Positions in Master Proteins" into plain vectors (no data.frame),
# so the per-row mapping loop does not build small tables.
pm_positions_parts <- function(x) {
  none <- list(accession = character(), start = integer(), end = integer(), parsed = TRUE)
  if (is.na(x) || !nzchar(trimws(x))) return(none)
  pattern <- "^(?:([^;[:space:]\\[\\]]+)[[:space:]]*)?\\[([0-9]+)-([0-9]+)\\]$"
  tokens <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  acc <- character(); st <- integer(); en <- integer(); parsed <- TRUE; accession <- NULL
  for (token in tokens[nzchar(tokens)]) {
    m <- regmatches(token, regexec(pattern, token, perl = TRUE))[[1]]
    if (!length(m)) {
      # Do not carry an accession across an unrecognized fragment.
      parsed <- FALSE; accession <- NULL; next
    }
    if (nzchar(m[2])) accession <- m[2]
    bounds <- as.numeric(m[3:4])
    if (is.null(accession) || any(!is.finite(bounds) | bounds > .Machine$integer.max)) {
      parsed <- FALSE; accession <- NULL; next
    }
    acc <- c(acc, accession); st <- c(st, as.integer(bounds[1])); en <- c(en, as.integer(bounds[2]))
  }
  list(accession = acc, start = st, end = en, parsed = parsed)
}
pm_positions <- function(x) {
  empty <- data.frame(accession = character(), start = integer(), end = integer())
  z <- pm_positions_parts(x)
  data <- if (length(z$accession)) data.frame(accession = z$accession, start = z$start, end = z$end) else empty
  list(data = data, parsed = z$parsed)
}
pm_occurrences <- function(sequence, protein) {
  # Fixed-string search advanced by one residue includes overlapping matches.
  out <- integer(); offset <- 0L
  while (offset <= nchar(protein) - nchar(sequence)) {
    p <- regexpr(sequence, substring(protein, offset + 1L), fixed = TRUE)[1]
    if (p < 0L) break
    start <- offset + p
    out <- c(out, start)
    offset <- start
  }
  out
}
# Complete analysis entrypoints reject malformed selected quantities. Low-level
# parsers still retain diagnostic statuses for callers explicitly inspecting them.
pm_validate_quantity_columns <- function(data, columns, label = "Peptide table") {
  for (column in unique(columns)) {
    if (!column %in% names(data)) stop("Unknown quantity column: ", column)
    original <- as.character(data[[column]])
    blank <- is.na(original) | trimws(original) %in% c("", "NA", "N/A", "n/a")
    invalid <- which(!blank & !is.finite(suppressWarnings(as.numeric(original))))
    if (length(invalid)) {
      shown <- head(invalid, 5L)
      examples <- paste(paste0("row ", shown + 1L, ": ", encodeString(original[shown], quote='"')), collapse="; ")
      stop(label, " quantity column '", column, "' contains ", length(invalid),
        " invalid numeric value(s) (rows include the header): ", examples,
        ". Use a decimal point (100.5), without thousands separators, or leave missing cells blank. No analysis results were produced.", call.=FALSE)
    }
  }
  invisible(TRUE)
}

pm_measurements <- function(x) {
  empty <- data.frame(peptide_id = character(), source_column = character(), abundance_type = character(),
                      measurement_label = character(), value = numeric(), original_value = character(), value_status = character())
  cols <- grep("^(Abundance:|Abundances \\(Normalized\\):|Abundances \\(Grouped\\):)", names(x), value = TRUE)
  rows <- lapply(cols, function(h) {
    original <- as.character(x[[h]])
    blank <- is.na(original) | trimws(original) %in% c("", "NA", "N/A", "n/a")
    value <- suppressWarnings(as.numeric(original))
    status <- ifelse(blank, "missing", ifelse(!is.finite(value), "invalid_numeric_value", "observed"))
    value[status != "observed"] <- NA_real_
    data.frame(peptide_id = sprintf("pep_%07d", seq_len(nrow(x))), source_column = h,
               abundance_type = if (startsWith(h, "Abundance:")) "abundance" else if (grepl("Normalized", h, fixed = TRUE)) "normalized" else "grouped",
               measurement_label = trimws(sub("^[^:]+:", "", h)), value = value,
               original_value = original, value_status = status, stringsAsFactors = FALSE)
  })
  pm_bind(rows, empty)
}
pm_empty_matches <- function() data.frame(peptide_id = character(), record_id = integer(), accession = character(),
  start = integer(), end = integer(), protein_length = integer(), left_flank = character(), right_flank = character(),
  is_master = logical(), reported_by_pd = logical(), pd_coordinate_agrees = logical(), stringsAsFactors = FALSE)
pm_map <- function(peptides, fasta, proteins = NULL, contaminant_ids = character(), contaminant_fasta = NULL,
                   scope = c("reported", "all"), sequence_column = "Annotated Sequence") {
  scope <- match.arg(scope)
  if (!sequence_column %in% names(peptides)) stop("Missing peptide sequence column: ", sequence_column)
  if (scope == "reported" && !any(c("Protein Accessions", "Master Protein Accessions", "Positions in Master Proteins") %in% names(peptides)))
    stop("Reported scope needs a PD accession column. Use scope='all' for FASTA-wide mapping.")
  idx <- list2env(split(seq_len(nrow(fasta)), fasta$accession), hash = TRUE, parent = emptyenv())
  lookup <- function(a) get0(a, idx, inherits = FALSE, ifnotfound = integer())
  contaminant_ids <- unique(contaminant_ids[nzchar(contaminant_ids)])
  if (!is.null(proteins) && all(c("Accession", "Contaminant") %in% names(proteins))) {
    flags <- pm_flag(proteins$Contaminant)
    contaminant_ids <- unique(c(contaminant_ids, as.character(proteins$Accession[which(flags %in% TRUE)])))
  }
  if (!is.null(contaminant_fasta)) contaminant_ids <- unique(c(contaminant_ids, contaminant_fasta$accession))
  contaminant_index <- which(fasta$accession %in% contaminant_ids)
  contaminant_sequences <- unique(c(fasta$sequence[contaminant_index], if (!is.null(contaminant_fasta)) contaminant_fasta$sequence else character()))
  flag_column_present <- "Contaminant" %in% names(peptides)
  flags <- if (flag_column_present) pm_flag(peptides$Contaminant) else rep(NA, nrow(peptides))
  seq_cache <- new.env(hash = TRUE, parent = emptyenv())
  contamination_cache <- new.env(hash = TRUE, parent = emptyenv())
  source_file <- attr(peptides, "source_file"); if (is.null(source_file)) source_file <- ""
  source_sheet <- attr(peptides, "source_sheet"); if (is.null(source_sheet)) source_sheet <- ""
  n <- nrow(peptides)
  if (!n) stop("Peptide table has no data rows.")
  # Per-row results accumulate in vectors/lists; tables are built once at the end
  # (same columns, types and order as the former per-row data.frame + rbind).
  column <- function(name) if (name %in% names(peptides)) as.character(peptides[[name]]) else rep(NA_character_, n)
  seq_v <- column(sequence_column); pos_v <- column("Positions in Master Proteins")
  master_v <- column("Master Protein Accessions"); accession_v <- column("Protein Accessions")
  ids <- sprintf("pep_%07d", seq_len(n))
  parsed_all <- pm_parse_sequence_all(seq_v)
  masters_all <- pm_split_all(master_v); accessions_all <- pm_split_all(accession_v)
  s_sequence <- s_left <- s_right <- rep(NA_character_, n)
  s_screening <- s_reason <- s_status <- character(n)
  s_count <- s_accessions <- integer(n); s_pd_valid <- logical(n)
  match_rows <- vector("list", n); checks <- vector("list", n); candidate_rows <- vector("list", n)
  raw_fields <- c("record_id","accession","start","end","protein_length","left_flank","right_flank")
  empty_raw <- as.list(pm_empty_matches()[raw_fields])
  for (i in seq_len(n)) {
    if (i %% 5000L == 0L) message("Mapping peptide row ", i, " / ", n)
    id <- ids[i]
    parsed <- list(sequence = parsed_all$sequence[i], left = parsed_all$left[i], right = parsed_all$right[i],
                   status = parsed_all$status[i]); sequence <- parsed$sequence
    pdp <- pm_positions_parts(pos_v[i])
    masters <- masters_all[[i]]
    reported <- unique(c(accessions_all[[i]], masters, pdp$accession))
    reported <- reported[nzchar(reported)]
    reason <- if (isTRUE(flags[i])) "pd_contaminant_flag" else ""
    listed <- reported %in% contaminant_ids
    if (!nzchar(reason) && any(listed)) reason <- if (all(listed)) "contaminant_accession" else "shared_contaminant_assignment"
    if (parsed$status == "parsed" && !nzchar(reason) && length(contaminant_sequences)) {
      if (!exists(sequence, contamination_cache, inherits = FALSE))
        assign(sequence, any(grepl(sequence, contaminant_sequences, fixed = TRUE)), contamination_cache)
      if (get(sequence, contamination_cache)) reason <- "contaminant_sequence_overlap"
    }
    screening <- if (nzchar(reason)) "excluded" else if (!is.na(flags[i])) "pd_flag_and_available_lists_checked" else if (length(contaminant_ids)) "available_lists_checked_pd_flag_unknown" else "not_available"
    nh <- 0L; hit_accessions <- character()
    # Contaminants stay in source/exclusion tables; no analysis matches are emitted.
    if (!nzchar(reason) && parsed$status == "parsed") {
      candidates <- if (scope == "all") seq_len(nrow(fasta)) else unique(unlist(lapply(reported, lookup), use.names = FALSE))
      if (is.null(candidates)) candidates <- integer()
      key <- paste(sequence, scope, paste(reported, collapse = ";"), sep = "|")
      if (!exists(key, seq_cache, inherits = FALSE)) {
        present <- candidates[grepl(sequence, fasta$sequence[candidates], fixed = TRUE)]
        raw_hits <- lapply(present, function(j) {
          starts <- pm_occurrences(sequence, fasta$sequence[j]); ends <- starts + nchar(sequence) - 1L
          list(record_id = rep(fasta$record_id[j], length(starts)), accession = rep(fasta$accession[j], length(starts)),
               start = starts, end = ends, protein_length = rep(fasta$length[j], length(starts)),
               left_flank = vapply(starts, function(s) if (s == 1L) "-" else substr(fasta$sequence[j], s-1L, s-1L), character(1)),
               right_flank = vapply(ends, function(e) if (e == fasta$length[j]) "-" else substr(fasta$sequence[j], e+1L, e+1L), character(1)))
        })
        assign(key, if (length(raw_hits)) setNames(lapply(raw_fields, function(f)
          unlist(lapply(raw_hits, `[[`, f), use.names = FALSE)), raw_fields) else empty_raw, seq_cache)
      }
      raw <- get(key, seq_cache)
      nh <- length(raw$start)
      if (nh) {
        agrees <- rep(NA, nh)
        for (h in seq_len(nh)) {
          same <- pdp$accession == raw$accession[h]
          if (any(same)) agrees[h] <- any(pdp$start[same] == raw$start[h] & pdp$end[same] == raw$end[h])
        }
        hit_accessions <- raw$accession
        match_rows[[i]] <- c(list(peptide_id = rep(id, nh)), raw,
          list(is_master = raw$accession %in% masters, reported_by_pd = raw$accession %in% reported, pd_coordinate_agrees = agrees))
      }
      if (length(reported)) {
        state <- vapply(reported, function(a) {
          ri <- lookup(a)
          if (!length(ri)) "accession_absent" else if (a %in% hit_accessions) "matched" else "sequence_absent"
        }, character(1), USE.NAMES = FALSE)
        candidate_rows[[i]] <- list(peptide_id = rep(id, length(reported)), accession = reported, status = state)
      }
      if (length(pdp$accession)) {
        state <- vapply(seq_along(pdp$accession), function(p) {
          pa <- pdp$accession[p]; ps <- pdp$start[p]; pe <- pdp$end[p]; ri <- lookup(pa)
          valid <- ps >= 1L && pe >= ps && pe - ps + 1L == nchar(sequence)
          if (!valid) "invalid_pd_interval" else if (!length(ri)) "accession_absent" else if (any(vapply(ri, function(j)
            pe <= fasta$length[j] && substr(fasta$sequence[j], ps, pe) == sequence, logical(1)))) "coordinate_match" else if (pa %in% hit_accessions) "different_position" else "sequence_absent"
        }, character(1))
        checks[[i]] <- list(peptide_id = rep(id, length(state)), accession = pdp$accession, pd_start = pdp$start,
                            pd_end = pdp$end, status = state)
      }
    }
    s_status[i] <- if (nzchar(reason)) "excluded_contaminant" else if (parsed$status != "parsed") parsed$status else if (!length(reported) && scope == "reported") "missing_reported_accessions" else if (!nh) "unmatched" else if (nh == 1L) "one_match_in_scope" else "multiple_matches_in_scope"
    s_sequence[i] <- sequence; s_left[i] <- parsed$left; s_right[i] <- parsed$right
    s_screening[i] <- screening; s_reason[i] <- reason
    s_count[i] <- nh; s_accessions[i] <- length(unique(hit_accessions)); s_pd_valid[i] <- pdp$parsed
  }
  summary <- data.frame(peptide_id = ids, source_file = rep(source_file, n), source_sheet = rep(source_sheet, n),
    source_row = seq_len(n) + 1L, annotated_sequence = seq_v, sequence = s_sequence,
    left_annotation = s_left, right_annotation = s_right, modifications = column("Modifications"),
    modifications_available = rep("Modifications" %in% names(peptides), n), contaminant_flag = flags,
    contaminant_screening = s_screening, exclusion_reason = s_reason, mapping_scope = rep(scope, n),
    mapping_status = s_status, match_count = s_count, matched_accession_count = s_accessions,
    pd_position_format_valid = s_pd_valid, stringsAsFactors = FALSE)
  # Combine per-row column lists into one table; `empty` supplies the zero-row schema.
  assemble <- function(rows, empty) {
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(empty)
    out <- lapply(names(empty), function(f) unlist(lapply(rows, `[[`, f), use.names = FALSE))
    names(out) <- names(empty)
    as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
  }
  matches <- assemble(match_rows, pm_empty_matches())
  pd_checks <- assemble(checks, data.frame(peptide_id = character(), accession = character(), pd_start = integer(), pd_end = integer(), status = character()))
  candidates <- assemble(candidate_rows, data.frame(peptide_id = character(), accession = character(), status = character()))
  measures <- pm_measurements(peptides)
  if (nrow(measures)) measures$included_in_analysis <- !measures$peptide_id %in% summary$peptide_id[summary$mapping_status == "excluded_contaminant"]
  list(peptides = summary, matches = matches, pd_checks = pd_checks, candidates = candidates, measurements = measures,
       excluded = summary[summary$mapping_status == "excluded_contaminant", , drop = FALSE],
       source_peptides = peptides, source_proteins = proteins, fasta = fasta,
       settings = list(scope = scope, sequence_column = sequence_column, exact_matching = TRUE,
                       il_equivalence = FALSE, coordinates = "one-based inclusive", contaminant_ids = contaminant_ids,
                       contamination_policy = "Contaminant inputs are optional. Exclude only on supplied flags or contaminant accession/sequence overlap; otherwise proceed without contaminant exclusions. Retain audit rows."))
}
