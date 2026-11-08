#!/usr/bin/env Rscript
# Run from any working directory: Rscript scripts/run.R --config config.json
args <- commandArgs(trailingOnly = TRUE)
# Windows R sessions started under a Unix-style locale can lose Unicode paths.
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  locale <- suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
  if (!nzchar(locale)) stop("A UTF-8 locale is required to read Unicode file paths.")
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)))
for (path in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
  source(path, encoding = "UTF-8")
if (length(args) != 2L || args[1] != "--config") stop("Usage: Rscript scripts/run.R --config config.json")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite first.")
config_path <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
config <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
pm_run(config, config_base = dirname(config_path), template = file.path(root,"inst","templates","report_template.html"))
