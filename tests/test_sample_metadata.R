args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]), winslash = "/")))
source(file.path(root, "R", "excel_locale.R"), encoding = "UTF-8")
source(file.path(root, "R", "sample_metadata.R"), encoding = "UTF-8")
stopifnot(requireNamespace("openxlsx", quietly = TRUE), requireNamespace("readxl", quietly = TRUE))
td <- tempfile("sample_metadata_"); dir.create(td)

# Multiple peptide rows retain one metadata row per actual measurement column.
raw <- c("Abundance: F1: Sample", "Abundance: F2: Sample", "Abundance: F3: Sample")
grouped <- c("Abundances (Grouped): F1", "Abundances (Grouped): F2")
measurements <- data.frame(peptide_id = rep(c("p1", "p2"), each = 5),
  source_column = rep(c(raw, grouped), 2), measurement_label = rep(c("F1: Sample", "F2: Sample", "F3: Sample", "F1", "F2"), 2),
  abundance_type = rep(c(rep("abundance", 3), rep("grouped", 2)), 2),
  value = c(0, 1, NA, 20, 30, 5, 6, 7, 20, 30), stringsAsFactors = FALSE)
proteins <- data.frame(Accession = c("P1", "P2"), check.names = FALSE)
proteins[[raw[1]]] <- c(10, NA); proteins[[raw[2]]] <- c(20, 0)
proteins[["Abundance: Protein F3"]] <- c(30, 40)
proteins[["Description"]] <- c("one", "two")
result <- list(measurements = measurements, source_proteins = proteins)
before <- serialize(result, NULL)
template <- pm_sample_template(result)
stopifnot(identical(names(template$samples), pm_sample_columns()), nrow(template$samples) == 3L,
  identical(template$samples[["Peptide column"]], raw), identical(template$samples[["Protein column"]], c(raw[1:2], "")),
  all(template$samples$Condition == ""), all(template$samples$Use == "Yes"), nrow(template$comparisons) == 1L,
  identical(pm_sample_template(result, "grouped")$samples[["Peptide column"]], grouped))
fails <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}
fails(pm_sample_template(result, "normalized"), "No measurement columns")
inconsistent <- result; inconsistent$measurements$measurement_label[6] <- "different label"
fails(pm_sample_template(inconsistent), "consistent source-column")

# A complete explicit mapping can use differently named protein measurements.
meta <- template
meta$samples$Sample <- c("00123", "\u00c7al\u0131\u015fma 2", "=1+1")
meta$samples$Condition <- c("Control", "Treatment", "Treatment")
meta$samples[["Protein column"]][3] <- "Abundance: Protein F3"
meta$comparisons <- data.frame(Name = "Treatment vs Control", Treatment = "Treatment", Control = "Control")
valid <- pm_validate_sample_metadata(meta, result, require_protein = TRUE)
stopifnot(identical(valid, meta), identical(before, serialize(result, NULL)))

# Import normalizes user input labels but never guesses source measurements.
normalized <- meta; normalized$samples$Use <- c(" yes ", "YES", "no")
normalized$samples$Condition[3] <- ""; normalized$samples[["Protein column"]][3] <- ""
normalized <- pm_validate_sample_metadata(normalized, result, TRUE)
stopifnot(identical(normalized$samples$Use, c("Yes", "Yes", "No")), normalized$samples$Condition[3] == "")
bad <- meta; bad$samples$Use[1] <- "Maybe"
fails(pm_validate_sample_metadata(bad, result), "Yes or No")
bad <- meta; bad$samples$Use <- "No"
fails(pm_validate_sample_metadata(bad, result), "at least one sample")
bad <- meta; bad$samples$Condition[1] <- ""
fails(pm_validate_sample_metadata(bad, result), "Sample and Condition")
bad <- meta; bad$samples$Sample[2] <- bad$samples$Sample[1]
fails(pm_validate_sample_metadata(bad, result), "identifiers must be unique")
bad <- meta; bad$samples[["Peptide column"]][2] <- raw[1]
fails(pm_validate_sample_metadata(bad, result), "Peptide column cannot be reused")
bad <- meta; bad$samples[["Peptide column"]][2] <- "Abundance: Missing"
fails(pm_validate_sample_metadata(bad, result), "Unknown Peptide column")
bad <- meta; bad$samples[["Peptide column"]][2] <- paste0(" ", raw[2])
fails(pm_validate_sample_metadata(bad, result), "Unknown Peptide column")
bad <- meta; bad$samples[["Peptide column"]][2] <- grouped[2]
fails(pm_validate_sample_metadata(bad, result), "one abundance type")
bad <- meta; bad$samples[["Protein column"]][2] <- raw[1]
fails(pm_validate_sample_metadata(bad, result, TRUE), "Protein column cannot be reused")
bad <- meta; bad$samples[["Protein column"]][2] <- "Abundance: Missing"
fails(pm_validate_sample_metadata(bad, result, TRUE), "Unknown Protein column")
bad <- meta; bad$samples[["Protein column"]][2] <- "Description"
fails(pm_validate_sample_metadata(bad, result, TRUE), "recognized abundance measurements")
bad <- meta; bad$samples[["Protein column"]][2] <- ""
fails(pm_validate_sample_metadata(bad, result, TRUE), "needs a Protein column")
stopifnot(nrow(pm_validate_sample_metadata(bad, result, FALSE)$samples) == 3L)
bad <- meta; bad$samples$Use[3] <- "No"; bad$samples$Condition[3] <- ""; bad$samples[["Protein column"]][3] <- "unknown unused header"
stopifnot(nrow(pm_validate_sample_metadata(bad, result, TRUE)$samples) == 3L)
bad <- meta; bad$comparisons <- template$comparisons
fails(pm_validate_sample_metadata(bad, result), "at least one comparison")
bad <- meta; bad$comparisons$Control <- "Unknown"
fails(pm_validate_sample_metadata(bad, result), "absent from included samples")
bad <- meta; bad$comparisons$Control <- "Treatment"
fails(pm_validate_sample_metadata(bad, result), "different conditions")
bad <- meta; bad$comparisons <- rbind(meta$comparisons, meta$comparisons)
fails(pm_validate_sample_metadata(bad, result), "names must be unique")
bad$comparisons$Name[2] <- "Duplicate contrast"
fails(pm_validate_sample_metadata(bad, result), "cannot be repeated")
bad <- meta; bad$samples$Extra <- "unexpected"
fails(pm_validate_sample_metadata(bad, result), "exactly these columns")

# Missing protein quantities do not block a peptide-only metadata file.
peptide_only <- result; peptide_only$source_proteins <- NULL
only <- pm_sample_template(peptide_only)
only$samples$Condition <- meta$samples$Condition; only$comparisons <- meta$comparisons
stopifnot(all(only$samples[["Protein column"]] == ""), nrow(pm_validate_sample_metadata(only, peptide_only)$samples) == 3L)
fails(pm_validate_sample_metadata(only, peptide_only, TRUE), "needs a Protein column")

# Saved XLSX roundtrip preserves literal identifiers, Unicode and exact headers.
path <- file.path(td, "metadata.xlsx")
pm_write_sample_metadata(meta, path)
roundtrip <- pm_read_sample_metadata(path, result, require_protein = TRUE)
stopifnot(identical(roundtrip, meta), identical(before, serialize(result, NULL)),
  identical(readxl::excel_sheets(path), c("Samples", "Comparisons", "Read me")))
fails(pm_write_sample_metadata(meta, path), "already exists")
unzip_dir <- file.path(td, "xlsx"); unzip(path, exdir = unzip_dir)
xml <- paste(readLines(file.path(unzip_dir, "xl", "worksheets", "sheet1.xml"), warn = FALSE), collapse = "")
stopifnot(grepl("dataValidation", xml, fixed = TRUE), grepl("Yes,No", xml, fixed = TRUE),
  grepl("frozen", xml, fixed = TRUE), !grepl("<f[ >]", xml))
blank_path <- file.path(td, "template.xlsx")
pm_write_sample_metadata(template, blank_path)
fails(pm_read_sample_metadata(blank_path, result), "Sample and Condition")
single <- meta; single$samples <- single$samples[1, , drop = FALSE]
fails(pm_validate_sample_metadata(single, result), "absent from included samples")
cat("Sample metadata tests passed.\n")
local({
  old<-options(warn=2);on.exit(options(old))
  known<-'one argument not used by format \'<x14:dataValidation type="list" allowBlank="%s" showInputMessage="%s" showErrorMessage="%s">\''
  withCallingHandlers(warning(known),warning=pm_openxlsx_list_warning)
  fails(withCallingHandlers(warning("Unrelated workbook warning"),warning=pm_openxlsx_list_warning),"Unrelated workbook warning")
})
cat("Only the known upstream list-validation warning is handled.\n")
