# Build the installable source archive from the repository root (a standard
# R package: DESCRIPTION, NAMESPACE, R/, man/, inst/templates/). Development
# files such as docs/, examples/, tests/ and tools/ are excluded by .Rbuildignore.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/build_package.R OUTPUT_DIRECTORY")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[1], mustWork = TRUE)
version <- read.dcf(file.path(root, "DESCRIPTION"))[1, "Version"]
archive <- file.path(out, paste0("atelente_", version, ".tar.gz"))
if (file.exists(archive)) stop("Archive exists. Choose a new output directory.")
assets <- file.path(root, "inst", "templates", c("report_template.html", "cleavage_report_template.html",
  "pair_summary.js", "report_payload.js", "report_workbook.js"))
if (!all(file.exists(assets))) stop("Missing report assets in inst/templates.")
setwd(out)
status <- system2(file.path(R.home("bin"), "R"), c("CMD", "build", "--no-build-vignettes", shQuote(root)))
if (status != 0L || !file.exists(archive)) stop("R CMD build failed.")
message("Built ", archive)
