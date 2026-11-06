# Repackaging a versioned compact report must not require expanded sample ratios.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash="/", mustWork=TRUE)))
if(.Platform$OS.type=="windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
for(f in list.files(file.path(root,"R"), "\\.R$", full.names=TRUE)) source(f, encoding="UTF-8")
fixture <- jsonlite::read_json(file.path(root,"tests/fixtures/export_bundle.json"), simplifyVector=TRUE)
td <- tempfile("atelente_compact_refresh_"); dir.create(td)
sha <- function(path) digest::digest(file=path, algo="sha256")
read_payload <- function(path) {
  html <- paste(readLines(path, encoding="UTF-8", warn=FALSE), collapse="\n")
  marker <- '<script type="application/json" id="data">'
  begin <- regexpr(marker, html, fixed=TRUE)[1]
  stopifnot(begin>0)
  tail <- substring(html, begin+nchar(marker), nchar(html))
  end <- regexpr('</script>', tail, fixed=TRUE)[1]
  stopifnot(end>0)
  pm_report_decode_json(substr(tail, 1, end-1L))
}
run_refresh <- function(input, name) {
  out <- file.path(td,name)
  log <- suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),
    shQuote(c(file.path(root,"tools/refresh_output_bundle.R"), input, out)), stdout=TRUE, stderr=TRUE))
  status <- attr(log,"status"); if(is.null(status)) status <- 0L
  if(status!=0L) stop(paste(log,collapse="\n"))
  out
}
read_sheets <- function(file) {
  names <- openxlsx::getSheetNames(file)
  setNames(lapply(names,function(name)as.data.frame(readxl::read_xlsx(file,
    sheet=name, trim_ws=FALSE, .name_repair="minimal"))),names)
}
manifest <- list(tool="Atelente",version="compact-refresh-fixture",
  settings=list(scope="reported"), metadata=list(note="Synthetic software fixture"))
baseline <- file.path(td,"source"); dir.create(baseline)
prepared <- pm_prepare_export_bundle(fixture,baseline,manifest,TRUE)
pm_write_report(prepared,file.path(root,"inst","templates","report_template.html"),file.path(baseline,"report.html"))
# A copied audit value is deliberately not used to rebuild the report: the
# compact payload itself carries the form measurements and prepared summaries.
sentinel <- file.path(baseline,"audit","source_sentinel.csv")
writeLines('"id","value"\n"001",1.23456789012345',sentinel,useBytes=TRUE)
manifest <- pm_finalize_audit(baseline,manifest)
jsonlite::write_json(manifest,file.path(baseline,"run_manifest.json"),auto_unbox=TRUE,
  pretty=TRUE,na="null",null="null",digits=NA)
source_files <- c(file.path(baseline,"report.html"),file.path(baseline,"atelente_results.xlsx"),
  file.path(baseline,"audit","atelente_tables.xlsx"),sentinel)
before <- setNames(vapply(source_files,sha,character(1)),source_files)
original <- read_payload(source_files[1])
unknown <- try(pm_report_decode_json('{"payload_schema":"atelente-report-data-999"}'),silent=TRUE)
stopifnot(inherits(unknown,"try-error"),grepl("Unsupported report data schema",as.character(unknown),fixed=TRUE))
stopifnot(identical(original$payload_schema,"atelente-report-data-2"),
  is.null(original$pair_view$ratios),is.null(original$pair_view$contrasts),
  is.data.frame(original$pair_view$prepared_summaries),
  nrow(original$pair_view$prepared_summaries)>0L,
  is.data.frame(original$pair_view$measurements),
  "__Results" %in% original$export_workbook$internal_sheets)
expected_public <- read_sheets(source_files[2])
expected_all <- read_sheets(source_files[3])
first <- run_refresh(source_files[1],"first_refresh")
second <- run_refresh(file.path(first,"report.html"),"second_refresh")
for(out in c(first,second)) {
  value <- read_payload(file.path(out,"report.html"))
  stopifnot(identical(value$payload_schema,original$payload_schema),
    is.null(value$pair_view$ratios),is.null(value$pair_view$contrasts),
    identical(value$pair_view,original$pair_view),
    identical(value$reference_access,original$reference_access),
    identical(read_sheets(file.path(out,"atelente_results.xlsx")),expected_public),
    identical(read_sheets(file.path(out,"audit","atelente_tables.xlsx")),expected_all),
    sha(file.path(out,"audit","source_sentinel.csv"))==before[[sentinel]])
}
stopifnot(identical(before,setNames(vapply(source_files,sha,character(1)),source_files)))
# Legacy reports remain readable and upgrade without needing a new analysis.
legacy <- file.path(td,"legacy"); dir.create(legacy)
json <- as.character(jsonlite::toJSON(fixture,dataframe="rows",auto_unbox=TRUE,
  digits=NA,na="null",null="null"))
writeLines(paste0('<script type="application/json" id="data">',json,'</script>'),
  file.path(legacy,"report.html"),useBytes=TRUE)
jsonlite::write_json(manifest,file.path(legacy,"run_manifest.json"),auto_unbox=TRUE,
  pretty=TRUE,na="null",null="null",digits=NA)
upgraded <- run_refresh(file.path(legacy,"report.html"),"legacy_refresh")
upgraded_payload <- read_payload(file.path(upgraded,"report.html"))
stopifnot(identical(upgraded_payload$payload_schema,"atelente-report-data-2"),
  is.null(upgraded_payload$pair_view$ratios),
  identical(read_sheets(file.path(upgraded,"atelente_results.xlsx")),expected_public),
  identical(read_sheets(file.path(upgraded,"audit","atelente_tables.xlsx")),expected_all))
node <- Sys.which("node")
if(nzchar(node)) {
  # Exercise the actual embedded application and its Current view callback,
  # then independently read the XLSX produced by the JavaScript exporter.
  log <- suppressWarnings(system2(node,shQuote(c(file.path(root,"tests/test_report_dom_smoke.cjs"),
    file.path(second,"report.html"))),stdout=TRUE,stderr=TRUE))
  status <- attr(log,"status"); if(is.null(status)) status <- 0L
  if(status!=0L) stop(paste(log,collapse="\n"))
  export_line <- log[grepl('^\\{"filtered_export":',log)]
  stopifnot(length(export_line)==1L)
  export_info <- jsonlite::fromJSON(export_line)
  view <- jsonlite::read_json(file.path(export_info$filtered_export,"expected.json"),simplifyVector=TRUE)
  output <- read_sheets(file.path(export_info$filtered_export,"filtered.xlsx"))
  index <- match(view$rows$key,expected_all$`__Results`$`Export key`)
  stopifnot(!anyNA(index),nrow(output$Results)>0L,nrow(output$Results)==length(index))
  for(column in intersect(c("Neo form ID","Spanning form ID","Neo log2FC","Spanning log2FC",
    "Ratio log2FC","P-value","BH-adjusted P","Test family","Test status"),names(output$Results))) {
    actual <- output$Results[[column]]; expected <- expected_all$`__Results`[[column]][index]
    # readxl infers all-blank filtered columns as logical even when the source
    # table contains numeric/text cells outside the exported subset.
    if(!identical(actual,expected) && !(length(actual)==length(expected)&&all(is.na(actual))&&all(is.na(expected))))
      stop("Current view column differs: ",column,"; ",paste(all.equal(actual,expected),collapse="; "))
  }
  stopifnot(identical(output$Results$Pattern,view$rows$pattern))
  cat("Compact application Current view Excel preserves saved effects/tests and click-time patterns.\n")
} else cat("Node unavailable: compact Current view application test not run.\n")
cat("Compact and legacy refresh: schema, repeated refresh, every workbook cell, audit bytes, source immutability and offline prepared summaries passed.\n")
