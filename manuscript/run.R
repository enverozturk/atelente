#!/usr/bin/env Rscript
# Run a prepared manuscript dataset through the installed Atelente package.
args <- commandArgs(trailingOnly = TRUE)
datasets <- c("01_GluC_PXD000503", "02_Kidney_MSV000094661", "03_Aortic_PXD073246")
if (length(args) != 1L || !args[1] %in% datasets)
  stop("Usage: Rscript manuscript/run.R DATASET; choose ", paste(datasets, collapse = ", "))
if (!requireNamespace("atelente", quietly = TRUE)) stop("Install Atelente using the repository README first.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
directory <- file.path(dirname(normalizePath(script, mustWork = TRUE)), args[1])
metadata <- jsonlite::fromJSON(file.path(directory, "metadata.json"))
config <- jsonlite::fromJSON(file.path(directory, "analysis.json"), simplifyVector = TRUE)
check_hash <- function(path, expected) {
  if (!identical(digest::digest(file = path, algo = "sha256"), expected))
    stop("Input checksum mismatch: ", basename(path))
}
peptides <- file.path(directory, "peptides.csv")
if (!file.exists(peptides)) {
  source <- gzfile(paste0(peptides, ".gz"), "rb")
  destination <- file(peptides, "wb")
  tryCatch({
    repeat {
      chunk <- readBin(source, what = "raw", n = 1048576L)
      if (!length(chunk)) break
      writeBin(chunk, destination)
    }
  }, finally = { close(source); close(destination) })
}
check_hash(peptides, metadata$input_sha256$peptides.csv)
# FASTA gzip is read directly by Atelente.
source <- gzfile(file.path(directory, config$fasta), "rb")
fasta <- tryCatch(readBin(source, what = "raw", n = 100000000L), finally = close(source))
if (!identical(digest::digest(fasta, algo = "sha256", serialize = FALSE), metadata$input_sha256$reference.fasta))
  stop("Input checksum mismatch: reference.fasta.gz")
if (isTRUE(config$pair_statistics)) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Install limma using the repository README first.")
  if (as.character(utils::packageVersion("limma")) != metadata$validated_limma_version)
    warning("The manuscript analyses used limma ", metadata$validated_limma_version,
            "; this run uses ", utils::packageVersion("limma"), ".")
}
atelente::pm_run(config, config_base = directory)
