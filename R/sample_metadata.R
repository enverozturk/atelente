# Explicit sample/condition metadata. This module never changes quantities.
pm_sample_columns <- function() c("Sample", "Condition", "Use", "Peptide column", "Protein column")
pm_comparison_columns <- function() c("Name", "Treatment", "Control")

pm_abundance_header_type <- function(headers) {
  out <- rep(NA_character_, length(headers))
  out[startsWith(headers, "Abundance:")] <- "abundance"
  out[startsWith(headers, "Abundances (Normalized):")] <- "normalized"
  out[startsWith(headers, "Abundances (Grouped):")] <- "grouped"
  out
}

pm_sample_template <- function(result, abundance_type = "abundance") {
  if (!is.character(abundance_type) || length(abundance_type) != 1L || is.na(abundance_type) || !nzchar(abundance_type))
    stop("Choose one abundance_type explicitly.")
  measures <- result$measurements
  needed <- c("source_column", "measurement_label", "abundance_type")
  if (is.null(measures) || !all(needed %in% names(measures)))
    stop("The result does not contain measurement-column metadata.")
  chosen <- measures[!is.na(measures$abundance_type) & measures$abundance_type == abundance_type,
                     needed, drop = FALSE]
  chosen <- unique(chosen)
  if (!nrow(chosen)) stop("No measurement columns found for abundance_type '", abundance_type, "'. Choose an available type explicitly.")
  if (anyNA(chosen$source_column) || any(!nzchar(chosen$source_column)) || anyDuplicated(chosen$source_column))
    stop("Measurement columns must have unique, consistent source-column metadata.")
  label <- as.character(chosen$measurement_label)
  label[is.na(label) | !nzchar(trimws(label))] <- chosen$source_column[is.na(label) | !nzchar(trimws(label))]
  protein_headers <- names(result$source_proteins)
  samples <- data.frame(Sample = label, Condition = rep("", nrow(chosen)), Use = rep("Yes", nrow(chosen)),
    "Peptide column" = as.character(chosen$source_column),
    "Protein column" = ifelse(chosen$source_column %in% protein_headers, chosen$source_column, ""),
    check.names = FALSE, stringsAsFactors = FALSE)
  comparisons <- data.frame(Name = "", Treatment = "", Control = "", stringsAsFactors = FALSE)
  list(samples = samples, comparisons = comparisons)
}

pm_metadata_table <- function(data, fields, label) {
  if (!is.data.frame(data)) stop(label, " must be a table.")
  if (anyDuplicated(names(data)) || !setequal(names(data), fields))
    stop(label, " must contain exactly these columns: ", paste(fields, collapse = ", "), ".")
  data <- data[, fields, drop = FALSE]
  for (name in fields) {
    if (is.list(data[[name]])) stop(label, " contains an unsupported column: ", name, ".")
    data[[name]] <- as.character(data[[name]])
    data[[name]][is.na(data[[name]])] <- ""
  }
  # Empty worksheet rows are layout, not samples or comparisons.
  data <- data[rowSums(do.call(cbind, lapply(data, function(x) nzchar(trimws(x))))) > 0L, , drop = FALSE]
  rownames(data) <- NULL
  data
}

pm_validate_sample_metadata <- function(metadata, result, require_protein = FALSE) {
  if (!is.logical(require_protein) || length(require_protein) != 1L || is.na(require_protein))
    stop("require_protein must be TRUE or FALSE.")
  samples <- pm_metadata_table(metadata$samples, pm_sample_columns(), "Samples")
  comparisons <- pm_metadata_table(metadata$comparisons, pm_comparison_columns(), "Comparisons")
  for (name in c("Sample", "Condition", "Use")) samples[[name]] <- trimws(samples[[name]])
  for (name in names(comparisons)) comparisons[[name]] <- trimws(comparisons[[name]])
  use <- tolower(samples$Use)
  if (any(!use %in% c("yes", "no"))) stop("Use must be Yes or No for every sample row.")
  samples$Use <- ifelse(use == "yes", "Yes", "No")
  selected <- samples$Use == "Yes"
  if (!any(selected)) stop("Include at least one sample with Use = Yes.")
  active <- samples[selected, , drop = FALSE]
  if (any(!nzchar(active$Sample)) || any(!nzchar(active$Condition)))
    stop("Every included sample needs Sample and Condition.")
  if (anyDuplicated(active$Sample)) stop("Included Sample identifiers must be unique.")
  if (any(!nzchar(active[["Peptide column"]]))) stop("Every included sample needs a Peptide column.")
  if (anyDuplicated(active[["Peptide column"]])) stop("A selected Peptide column cannot be reused across samples.")
  measures <- result$measurements
  needed <- c("source_column", "abundance_type")
  if (is.null(measures) || !all(needed %in% names(measures))) stop("The result does not contain measurement-column metadata.")
  lookup <- unique(measures[, needed, drop = FALSE])
  if (anyDuplicated(lookup$source_column)) stop("Measurement source columns have inconsistent abundance types.")
  index <- match(active[["Peptide column"]], lookup$source_column)
  if (anyNA(index)) stop("Unknown Peptide column: ", paste(active[["Peptide column"]][is.na(index)], collapse = ", "), ".")
  types <- unique(as.character(lookup$abundance_type[index]))
  if (length(types) != 1L || anyNA(types)) stop("Included peptide columns must use one abundance type; do not mix raw, normalized and grouped measurements.")

  # Explicit protein columns are checked even when correction is not requested.
  # Blank columns can remain in a peptide-only metadata workbook.
  protein <- active[["Protein column"]]
  supplied <- nzchar(protein)
  if (require_protein && any(!supplied)) stop("Every included sample needs a Protein column for protein correction.")
  if (anyDuplicated(protein[supplied])) stop("A selected Protein column cannot be reused across samples.")
  if (any(supplied)) {
    protein_headers <- names(result$source_proteins)
    unknown <- supplied & !protein %in% protein_headers
    if (any(unknown)) stop("Unknown Protein column: ", paste(protein[unknown], collapse = ", "), ".")
    protein_types <- pm_abundance_header_type(protein[supplied])
    if (anyNA(protein_types)) stop("Protein columns must refer to recognized abundance measurements.")
    if (length(unique(protein_types)) != 1L) stop("Included protein columns must use one abundance type.")
  }
  if (!nrow(comparisons)) stop("Add at least one comparison with Name, Treatment and Control.")
  if (any(!nzchar(comparisons$Name)) || any(!nzchar(comparisons$Treatment)) || any(!nzchar(comparisons$Control)))
    stop("Every comparison needs Name, Treatment and Control.")
  if (anyDuplicated(comparisons$Name)) stop("Comparison names must be unique.")
  if (any(comparisons$Treatment == comparisons$Control)) stop("Treatment and Control must be different conditions.")
  conditions <- unique(active$Condition)
  unknown <- setdiff(unique(c(comparisons$Treatment, comparisons$Control)), conditions)
  if (length(unknown)) stop("Comparison conditions are absent from included samples: ", paste(unknown, collapse = ", "), ".")
  if (anyDuplicated(comparisons[, c("Treatment", "Control"), drop = FALSE]))
    stop("The same Treatment/Control comparison cannot be repeated.")
  list(samples = samples, comparisons = comparisons)
}

# openxlsx 4.2.5.2 passes an unused argument to sprintf for list validation.
# The resulting XML is unchanged. Muffle this exact upstream message only;
# every other warning remains visible (including under options(warn=2)).
pm_openxlsx_list_warning <- function(w) {
  expected <- 'one argument not used by format \'<x14:dataValidation type="list" allowBlank="%s" showInputMessage="%s" showErrorMessage="%s">\''
  if (identical(conditionMessage(w), expected)) invokeRestart("muffleWarning")
}

pm_write_sample_metadata <- function(metadata, path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Install openxlsx to write sample metadata.")
  if (file.exists(path)) stop("Sample metadata already exists; choose a new path to preserve edits.")
  # Incomplete templates are intentionally writable; validation occurs on import.
  tables <- list(Samples = pm_metadata_table(metadata$samples, pm_sample_columns(), "Samples"),
                 Comparisons = pm_metadata_table(metadata$comparisons, pm_comparison_columns(), "Comparisons"))
  pm_with_excel_utf8({
    wb <- openxlsx::createWorkbook(creator = "Atelente", title = "Sample and condition mapping")
    openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 11, fontColour = "#172A39")
    header <- openxlsx::createStyle(fontColour = "#FFFFFF", fgFill = "#263F50", textDecoration = "bold", wrapText = TRUE, valign = "center")
    input <- openxlsx::createStyle(fontColour = "#154C9A", fgFill = "#EAF2FB", numFmt = "@", valign = "center")
    source <- openxlsx::createStyle(fontColour = "#526674", fgFill = "#F1F4F5", numFmt = "@", wrapText = TRUE, valign = "center")
    for (name in names(tables)) {
      data <- tables[[name]]
      openxlsx::addWorksheet(wb, name, gridLines = FALSE)
      openxlsx::writeData(wb, name, data, rowNames = FALSE, colNames = TRUE, headerStyle = header, withFilter = nrow(data) > 0L, keepNA = FALSE)
      openxlsx::freezePane(wb, name, firstRow = TRUE, firstCol = TRUE)
      openxlsx::setRowHeights(wb, name, 1, 28)
      last <- max(2L, nrow(data) + 1L)
      openxlsx::addStyle(wb, name, input, rows = 2:last, cols = 1:3, gridExpand = TRUE)
      openxlsx::setColWidths(wb, name, 1:3, widths = if (name == "Samples") c(25, 25, 12) else c(32, 25, 25))
      openxlsx::setRowHeights(wb, name, 2:last, 30)
      if (name == "Samples") {
        openxlsx::addStyle(wb, name, source, rows = 2:last, cols = 4:5, gridExpand = TRUE)
        openxlsx::setColWidths(wb, name, 4:5, widths = 48)
        withCallingHandlers(openxlsx::dataValidation(wb, name, cols = 3, rows = 2:last, type = "list", value = '"Yes,No"', allowBlank = FALSE,
          showInputMsg = FALSE, showErrorMsg = TRUE), warning=pm_openxlsx_list_warning)
      }
    }
    openxlsx::addWorksheet(wb, "Read me", gridLines = FALSE)
    notes <- c("Sample and condition mapping",
      "Samples: edit the blue Sample, Condition and Use cells. Use No for channels outside this comparison.",
      "Gray measurement columns preserve source headers. Change a mapping only when selecting the corresponding measurement explicitly.",
      "Comparisons: enter a name and the exact Treatment and Control condition labels. Positive log2 fold change means higher in Treatment.",
      "Each included row represents one measurement. Raw, normalized and grouped columns must not be mixed or counted as separate repeats.",
      "Protein column may remain blank for peptide-only use. Protein correction requires a matching abundance column for every included sample.")
    openxlsx::writeData(wb, "Read me", data.frame(Notes = notes), colNames = FALSE)
    openxlsx::setColWidths(wb, "Read me", 1, 100)
    openxlsx::setRowHeights(wb, "Read me", 1, 28)
    openxlsx::setRowHeights(wb, "Read me", 2:length(notes), 36)
    openxlsx::addStyle(wb, "Read me", header, rows = 1, cols = 1)
    openxlsx::addStyle(wb, "Read me", openxlsx::createStyle(wrapText = TRUE, valign = "center"), rows = 2:length(notes), cols = 1, gridExpand = TRUE)
    openxlsx::saveWorkbook(wb, path, overwrite = FALSE)
  })
  invisible(path)
}

pm_read_sample_metadata <- function(path, result, require_protein = FALSE) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read sample metadata.")
  sheets <- readxl::excel_sheets(path)
  if (!all(c("Samples", "Comparisons") %in% sheets)) stop("Sample metadata requires Samples and Comparisons worksheets.")
  read <- function(sheet) as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text", trim_ws = FALSE, .name_repair = "minimal"), stringsAsFactors = FALSE, check.names = FALSE)
  pm_validate_sample_metadata(list(samples = read("Samples"), comparisons = read("Comparisons")), result, require_protein)
}
