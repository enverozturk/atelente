# User-facing mapping workflow. Low-level mapping/export functions remain usable
# independently; the complete workflow requests reference annotations by default.
pm_run <- function(config, config_base = getwd(), template = NULL) {
  if (!is.list(config) || is.null(names(config)) || anyDuplicated(names(config)))
    stop("config must be a named list with unique fields.")
  config_base <- normalizePath(config_base, winslash = "/", mustWork = TRUE)
  scalar_flag <- function(x, name, default) {
    if (is.null(x)) return(default)
    if (!is.logical(x) || length(x) != 1L || is.na(x))
      stop(name, " must be TRUE or FALSE.")
    x
  }
  annotations <- scalar_flag(config$annotations, "annotations", TRUE)
  excel <- scalar_flag(config$excel, "excel", TRUE)
  pair_statistics <- scalar_flag(config$pair_statistics, "pair_statistics", FALSE)
  assignment_mode <- if (is.null(config$assignment_mode))
    if (is.null(config$column_mapping)) "unique" else "representative" else config$assignment_mode
  if (!is.character(assignment_mode) || length(assignment_mode) != 1L ||
      is.na(assignment_mode) || !assignment_mode %in% c("unique", "representative"))
    stop("assignment_mode must be unique or representative.")
  if (!is.null(config$sample_design) && !is.null(config$neo_spanning_samples))
    stop("Choose inline sample_design or neo_spanning_samples workbook, not both.")
  if (!is.null(config$sample_design) && is.null(config$column_mapping))
    stop("Inline sample_design requires explicit column_mapping.")
  if (!is.null(config$column_mapping) && !is.null(config$neo_spanning_samples))
    stop("Column-mapped input requires inline sample_design rather than a source-header workbook.")
  if (excel && !requireNamespace("openxlsx", quietly = TRUE))
    stop('Install openxlsx for Excel export, or set "excel": false in the config.')
  required <- c("peptides", "fasta", "output_dir")
  if (!all(required %in% names(config)))
    stop("Config requires peptides, fasta and output_dir.")
  resolve <- function(path) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
      stop("Configuration paths must be non-empty character scalars.")
    path <- path.expand(path)
    if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) file.path(config_base, path) else path
  }
  for (name in intersect(c(required, "proteins", "contaminant_ids", "contaminant_fasta",
                          "neo_spanning_samples", "cache_dir"), names(config)))
    config[[name]] <- resolve(config[[name]])
  for (name in intersect(c("merops", "uniprot"), names(config))) {
    spec <- config[[name]]
    if (!is.list(spec) || !all(c("snapshot", "provenance") %in% names(spec)))
      stop(name, " requires snapshot and provenance paths.")
    config[[name]] <- list(snapshot = resolve(spec$snapshot), provenance = resolve(spec$provenance))
  }
  if (is.null(template))
    template <- system.file("templates", "report_template.html", package = "atelente", mustWork = TRUE)
  if (!file.exists(template)) stop("Report template does not exist.")
  if (dir.exists(config$output_dir) && length(list.files(config$output_dir, all.files = TRUE, no.. = TRUE)))
    stop("Choose an empty output directory.")

  message("Reading peptide table and FASTA...")
  peptides <- pm_read_table(config$peptides, config$peptide_sheet)
  prepared <- NULL
  if (!is.null(config$column_mapping)) {
    prepared <- pm_prepare_column_input(peptides, config$column_mapping,
      config$sample_design, assignment_mode, pair_statistics)
    source_samples <- prepared$mapping_audit$source_samples
    pm_validate_quantity_columns(peptides,
      source_samples[["Peptide column"]][source_samples$Use == "Yes"])
    peptides <- prepared$peptides
    if (!is.null(config$sequence_column) && !identical(config$sequence_column, "Sequence"))
      stop("sequence_column conflicts with canonical Sequence from column_mapping.")
    config$sequence_column <- "Sequence"
  }
  representative_column <- if (!is.null(prepared)) "Master Protein Accessions" else config$representative_column
  if (assignment_mode == "representative" &&
      (!is.character(representative_column) || length(representative_column) != 1L ||
       is.na(representative_column) || !representative_column %in% names(peptides)))
    stop("Representative assignment requires an explicitly selected representative_column.")
  proteins <- if (!is.null(config$proteins)) pm_read_table(config$proteins, config$protein_sheet) else NULL
  quantity_headers <- function(x) grep("^(Abundance:|Abundances \\(Normalized\\):|Abundances \\(Grouped\\):)", names(x), value=TRUE)
  if (is.null(prepared)) pm_validate_quantity_columns(peptides, quantity_headers(peptides))
  if (!is.null(proteins)) pm_validate_quantity_columns(proteins, quantity_headers(proteins), "Protein table")
  fasta <- pm_read_fasta(config$fasta)
  contaminants <- if (!is.null(config$contaminant_ids))
    trimws(readLines(config$contaminant_ids, encoding = "UTF-8", warn = FALSE)) else character()
  contaminants <- contaminants[nzchar(contaminants) & !startsWith(contaminants, "#")]
  contaminant_fasta <- if (!is.null(config$contaminant_fasta)) pm_read_fasta(config$contaminant_fasta) else NULL
  scope <- if (is.null(config$scope)) "reported" else config$scope
  sequence_column <- if (is.null(config$sequence_column)) "Annotated Sequence" else config$sequence_column
  result <- pm_map(peptides, fasta, proteins, contaminants, contaminant_fasta, scope, sequence_column)
  if (!is.null(config$assignment_mode) || !is.null(prepared) || pair_statistics)
    result$settings$assignment_mode <- assignment_mode
  if (!is.null(prepared)) {
    result$source_input <- prepared$original
    result$column_mapping_audit <- prepared$mapping_audit
  }
  if (!is.null(config$cleavage)) {
    message("Linking candidate cleavage evidence...")
    result <- pm_build_cleavage(result, config$cleavage, digest::digest(file=config$fasta, algo="sha256"))
  }
  input_names <- intersect(c("peptides", "proteins", "fasta", "contaminant_ids", "contaminant_fasta"), names(config))
  design <- NULL
  experimental_units <- config$experimental_units
  if (!is.null(prepared) && nrow(prepared$samples)) {
    design <- prepared$design
    experimental_units <- prepared$experimental_units
    if (!is.null(config$experimental_units))
      stop("Declare experimental units in inline sample_design only, not twice.")
    # The adapter has validated columns, sample identities and optional contrasts.
    # A contrast-free run keeps explicit sample names without inventing conditions.
    if (!is.null(config$sample_abundance_type) && !identical(config$sample_abundance_type, "abundance"))
      stop("Column-mapped individual quantities use canonical abundance columns.")
    config$sample_abundance_type <- "abundance"
  }
  if (!is.null(config$neo_spanning_samples)) {
    if (is.null(result$cleavage)) stop("neo_spanning_samples requires a cleavage profile.")
    design <- pm_read_sample_metadata(config$neo_spanning_samples, result)
    active <- design$samples[design$samples$Use == "Yes", , drop = FALSE]
    type <- unique(pm_abundance_header_type(active[["Peptide column"]]))
    if (length(type) != 1L || !type %in% c("abundance", "normalized"))
      stop("Neo/spanning contrasts require individual sample columns of one abundance type.")
    if (!is.null(config$sample_abundance_type) && !identical(type, config$sample_abundance_type))
      stop("sample_abundance_type conflicts with the neo/spanning sample workbook.")
    config$sample_abundance_type <- type
    input_names <- c(input_names, "neo_spanning_samples")
  }
  if (pair_statistics) {
    if (is.null(result$cleavage) || is.null(design) || !nrow(design$comparisons))
      stop("Pair tests require a cleavage profile and an explicit sample comparison.")
    if (is.null(experimental_units)) stop("Pair tests require explicit independent experimental units.")
  }
  descriptive_assignment <- assignment_mode == "representative" &&
    !is.null(result$cleavage) && !is.null(design) && nrow(design$comparisons) > 0L &&
    nrow(result$measurements) > 0L
  if (pair_statistics || descriptive_assignment) {
    if (!exists("pm_pair_workflow", mode = "function")) stop("Pair workflow is unavailable.")
    # Assignment gates apply to descriptive comparisons too. Inference remains
    # explicitly optional and alone requires independent experimental units.
    result$pair_statistics <- pm_pair_workflow(result = result, design = design,
      test = pair_statistics, experimental_units = if (pair_statistics) experimental_units else NULL,
      assignment_mode = assignment_mode,
      representative_column = representative_column, abundance_type = config$sample_abundance_type)
  }

  result <- pm_annotate_references(result, annotations = annotations,
    merops = config$merops, uniprot = config$uniprot, cache_dir = config$cache_dir)
  manifest <- pm_export(result, config$output_dir, template, config[input_names],
    metadata = if (is.null(config$metadata)) list() else config$metadata, excel = excel,
    sample_abundance_type = config$sample_abundance_type, neo_design = design,
    representative_selection_k = config$representative_selection_k, compact_bundle = TRUE)
  message("Done. Report: ", file.path(config$output_dir, "report.html"))
  invisible(manifest)
}
