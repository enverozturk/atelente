# Optional argument verifies the same contract through an installed package.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash="/", mustWork=TRUE)))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) <= 1L)
if (length(args)) {
  .libPaths(c(normalizePath(args[1], mustWork=TRUE), .libPaths()))
  api <- loadNamespace("atelente", lib.loc=.libPaths()[1])
  template <- NULL
} else {
  api <- new.env(parent=globalenv())
  for (f in list.files(file.path(root,"R"), "\\.R$", full.names=TRUE))
    sys.source(f, envir=api)
  template <- file.path(root,"inst","templates","report_template.html")
}
td <- tempfile("mapping_only_export_")
config <- api$pm_wizard_demo(td)
config$excel <- TRUE
input_path <- file.path(td,"peptides.csv")
source <- api$pm_read_table(input_path)
source[["Signal 1"]][1:4] <- c("0", NA_character_, "1.25e-12", "23.456789012345")
api$pm_write_csv(source,input_path)
input_hash <- digest::digest(file=input_path,algo="sha256")

read_payload <- function(path) {
  html <- paste(readLines(path,encoding="UTF-8",warn=FALSE),collapse="\n")
  json <- strsplit(html,'<script type="application/json" id="data">',fixed=TRUE)[[1]][2]
  json <- strsplit(json,"</script>",fixed=TRUE)[[1]][1]
  api$pm_report_decode_json(json)
}
verify <- function(name, cleavage=FALSE, compare=TRUE) {
  run <- config
  if (!cleavage) run$cleavage <- NULL
  if (!compare) {
    run$sample_design$comparisons <- NULL
    run$sample_design$samples$Condition <- ""
  }
  run$output_dir <- name
  api$pm_run(run,config_base=td,template=template)
  output <- file.path(td,name)
  workbook <- file.path(output,"atelente_results.xlsx")
  peptides <- as.data.frame(readxl::read_excel(workbook,sheet="Peptides"),check.names=FALSE)
  settings <- as.data.frame(readxl::read_excel(workbook,sheet="Settings"),check.names=FALSE)
  samples <- settings[settings$Section=="Samples",,drop=FALSE]
  quantities <- paste0("Quantity: ",config$sample_design$samples$Sample)
  audit <- api$pm_read_table(file.path(output,"audit","measurements.csv"))
  payload <- read_payload(file.path(output,"report.html"))
  stopifnot(nrow(peptides)==6L,nrow(audit)==36L,nrow(samples)==6L,
    all(quantities %in% names(peptides)),
    identical(samples$Sample,config$sample_design$samples$Sample),
    identical(samples$`Source column`,paste0("Signal ",1:6)),
    all(samples$`Input scale`=="linear"),all(samples$`Export scale`=="linear"),
    identical(digest::digest(file=input_path,algo="sha256"),input_hash))
  for (i in 1:6) {
    expected <- as.numeric(source[[paste0("Signal ",i)]][peptides$`Source row`-1L])
    actual <- peptides[[quantities[i]]]
    stopifnot(is.numeric(actual),identical(is.na(actual),is.na(expected)),
      isTRUE(all.equal(actual,expected,tolerance=1e-13)))
  }
  source_rows <- match(2:4,peptides$`Source row`)
  stopifnot(!anyNA(source_rows),peptides[[quantities[1]]][source_rows[1]]==0,
    is.na(peptides[[quantities[1]]][source_rows[2]]),
    peptides[[quantities[1]]][source_rows[3]]>0)
  if (cleavage) {
    stopifnot(!is.null(payload$cleavage),!is.null(payload$pair_view),
      is.null(payload$sample_quantities),nrow(payload$pair_view$measurements)==36L)
  } else {
    stopifnot(is.null(payload$cleavage),is.null(payload$pair_view),
      nrow(payload$sample_quantities$measurements)==36L,
      !any(c("ratios","contrasts","selection","selection_k") %in% names(payload$sample_quantities)),
      !any(settings$Setting %in% c("Representative selection","Frozen selection k","Default display coverage k")),
      !any(grepl("neo_spanning_(ratios|contrasts)",list.files(file.path(output,"audit")))))
  }
  if (!compare) stopifnot(!any(grepl("log2FC",names(peptides))),
    all(is.na(samples$Condition)|samples$Condition==""))
}
verify("mapping_with_comparison")
verify("mapping_without_comparison",compare=FALSE)
verify("cleavage_with_comparison",cleavage=TRUE)
cat("Mapping-only compact quantities, missing/zero values and sample provenance passed; cleavage exports retained.\n")
