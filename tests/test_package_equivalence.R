# Install the built archive in a separate library before running this test.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tests/test_package_equivalence.R PACKAGE_LIBRARY")
.libPaths(c(normalizePath(args[1], mustWork = TRUE), .libPaths()))
ns <- loadNamespace("atelente", lib.loc = .libPaths()[1])
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
legacy <- new.env(parent = globalenv())
for (path in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) {
  sys.source(path, envir = legacy, keep.source = FALSE)
}
td <- tempfile("atelente-equivalence-")
dir.create(td)
ref <- file.path(td, "reference.fasta")
writeLines(c(">sp|P1|SYNTHETIC GN=TEST", "MACDEFGHIK"), ref)
p <- data.frame("Annotated Sequence" = c("MACD", "EFGHIK", "CDEFGH", "EFGHIK"),
                "Master Protein Accessions" = "P1",
                Modifications = c("", "", "", "1xOxidation"),
                "Abundance: F1: Sample" = c(90, 100, 50, NA),
                "Abundance: F2: Sample" = c(70, 80, 0, 40), check.names = FALSE)
profile <- list(enzyme = "trypsin", lysine_blocked = FALSE)
sha <- digest::digest(file = ref, algo = "sha256")
run_core <- function(api) {
  x <- api$pm_map(p, api$pm_read_fasta(ref))
  api$pm_build_cleavage(x, profile, sha)
}
before <- run_core(legacy)
after <- run_core(ns)
stopifnot(identical(before, after))
ratios <- ns$pm_neo_spanning_ratios(after$cleavage$evidence, after$measurements)
stopifnot(nrow(ratios) > 0L, any(is.finite(ratios$ratio)), anyNA(ratios$ratio),
          identical(ratios, legacy$pm_neo_spanning_ratios(before$cleavage$evidence, before$measurements)))
for (view in c("mapping", "cleavage")) {
  x <- if (view == "mapping") ns$pm_map(p, ns$pm_read_fasta(ref)) else after
  old <- file.path(td, paste0(view, "-source"))
  new <- file.path(td, paste0(view, "-installed"))
  inputs <- list(fasta = ref)
  legacy$pm_export(x, old, file.path(root,"inst","templates","report_template.html"), inputs)
  ns$pm_export(x, new, inputs = inputs)
  stopifnot(identical(list.files(old), list.files(new)))
  files <- list.files(old, pattern = "\\.(csv|html)$")
  stopifnot(identical(unname(tools::md5sum(file.path(old, files))),
                      unname(tools::md5sum(file.path(new, files)))))
  manifest <- function(dir) {
    m <- jsonlite::read_json(file.path(dir, "run_manifest.json"))
    m$created_utc <- NULL
    m
  }
  stopifnot(identical(manifest(old), manifest(new)))
}
if (requireNamespace("limma", quietly = TRUE)) {
  samples <- data.frame(sample = paste0("s", 1:8), condition = rep(c("WT", "KO"), each = 4),
                        experimental_unit = paste0("u", 1:8))
  response <- matrix(sin(seq_len(64L)), ncol = 8L,
                     dimnames = list(paste0("p", 1:8), samples$sample))
  response[1, 1] <- NA_real_
  family <- data.frame(measurement_pair_id = rownames(response), family = "N")
  stopifnot(identical(ns$pm_test_pair_matrix(response, samples, family, "WT", "KO"),
                      legacy$pm_test_pair_matrix(response, samples, family, "WT", "KO")))
}
verify_encoding <- function() {
  original <- Sys.getlocale("LC_CTYPE")
  on.exit(Sys.setlocale("LC_CTYPE", original))
  stopifnot(Sys.setlocale("LC_CTYPE", "C") == "C")
  path <- file.path(td, "utf8.csv")
  value <- intToUtf8(c(231, 305, 287))
  tab <- setNames(data.frame(note = c(value, "next row")), value)
  ns$pm_write_csv(tab, path)
  bytes <- readBin(path, "raw", file.info(path)$size)
  restored <- ns$pm_read_table(path)
  stopifnot(identical(names(restored), names(tab)), identical(restored[[1]], tab[[1]]),
            identical(Sys.getlocale("LC_CTYPE"), "C"),
            identical(bytes[1:3], as.raw(c(239, 187, 191))),
            length(grepRaw(charToRaw(value), bytes, fixed = TRUE)) > 0L)
}
verify_encoding()
# Installed R generation must also embed the lossless decoder and compressed data.
large <- list(cleavage=list(), text=paste(rep("UTF-8 \u03b1 test ", 100000),collapse=""))
source_report <- file.path(td,"compressed-source.html")
installed_report <- file.path(td,"compressed-installed.html")
legacy$pm_write_report(large,file.path(root,"inst","templates","cleavage_report_template.html"),source_report)
ns$pm_write_report(large,system.file("templates","cleavage_report_template.html",package="atelente"),installed_report)
stopifnot(identical(unname(tools::md5sum(source_report)),unname(tools::md5sum(installed_report))))
html <- paste(readLines(installed_report,encoding="UTF-8",warn=FALSE),collapse="\n")
stopifnot(grepl("atelente-report-envelope-1",html,fixed=TRUE),!grepl("__ATELENTE_REPORT_PAYLOAD_JS__",html,fixed=TRUE))
embedded <- sub('.*<script type="application/json" id="data">',"",html)
embedded <- substr(embedded,1,regexpr("</script>",embedded,fixed=TRUE)[1]-1L)
decoded <- ns$pm_report_decode_json(embedded,simplifyVector=FALSE)
stopifnot(identical(decoded$payload_schema,"atelente-report-data-2"))
decoded$payload_schema <- NULL
stopifnot(identical(decoded,large))
cat("Installed/source equivalence and UTF-8 CSV checks passed.\n")
