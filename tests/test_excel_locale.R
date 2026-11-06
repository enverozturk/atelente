# Excel writing must preserve non-ASCII text when R starts in a non-UTF-8
# locale (LC_ALL=C on Linux, Docker, CI). Unicode literals use \u escapes so
# this test does not depend on how its own source file is decoded.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
source(file.path(root,"R","excel_locale.R"),encoding="UTF-8")
source(file.path(root,"R","excel.R"),encoding="UTF-8")
if(!requireNamespace("openxlsx",quietly=TRUE) || !requireNamespace("readxl",quietly=TRUE)) {
  cat("openxlsx/readxl unavailable; Excel locale test skipped.\n"); quit(save="no")
}
td <- tempfile("excel-locale-"); dir.create(td)

# 1. Restoration and refusal, checked in-process from a C character type.
original <- Sys.getlocale("LC_CTYPE"); numeric <- Sys.getlocale("LC_NUMERIC")
if(!nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE","C")))) stop("Cannot set LC_CTYPE=C for the test.")
on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE",original)),add=TRUE)
stopifnot(!isTRUE(l10n_info()[["UTF-8"]]))
seen <- pm_with_excel_utf8(c(ctype=Sys.getlocale("LC_CTYPE"),utf8=isTRUE(l10n_info()[["UTF-8"]])))
stopifnot(seen[["utf8"]]=="TRUE", Sys.getlocale("LC_CTYPE")=="C", Sys.getlocale("LC_NUMERIC")==numeric)
# Restored after an error raised inside the UTF-8 block.
stopifnot(inherits(try(pm_with_excel_utf8(stop("deliberate failure")),silent=TRUE),"try-error"),
          Sys.getlocale("LC_CTYPE")=="C", Sys.getlocale("LC_NUMERIC")==numeric)
# Without a usable UTF-8 locale: clear error, nothing written, locale unchanged.
never <- file.path(td,"never.xlsx")
refused <- try(pm_with_excel_utf8(file.create(never),candidates="xx_NOT_A_LOCALE.UTF-8"),silent=TRUE)
stopifnot(inherits(refused,"try-error"), grepl("UTF-8 locale",refused,fixed=TRUE),
          !file.exists(never), Sys.getlocale("LC_CTYPE")=="C")
# An already UTF-8 session is left untouched (NULL means "not changed").
if(nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE","C.UTF-8"))) && isTRUE(l10n_info()[["UTF-8"]])) {
  stopifnot(is.null(pm_excel_utf8_ctype()), Sys.getlocale("LC_CTYPE")=="C.UTF-8")
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE","C")))
}

# 2. Full workbook writes in separate R processes started under C and UTF-8.
child <- file.path(td,"child.R")
writeLines(c(
  sprintf('source(%s,encoding="UTF-8");source(%s,encoding="UTF-8")',
          deparse(file.path(root,"R","excel_locale.R")),deparse(file.path(root,"R","excel.R"))),
  'args <- commandArgs(TRUE)',
  'stopifnot(nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE",args[2]))))',
  'before <- Sys.getlocale("LC_CTYPE")',
  'label <- c("\\u0130la\\u00e7 \\u03b1 & <1>","\\u96f6 \\U0001F9EA","plain")',
  'x <- data.frame(a=label,b=c(1.23456789012345e-120,NA,0),c=c(TRUE,NA,FALSE),stringsAsFactors=FALSE)',
  'names(x) <- c(label[1],label[2],"Flag")',
  'pm_write_excel_tables(list(Results=x,Settings=data.frame(Key=label[2],Value=label[1])),args[1])',
  'stopifnot(Sys.getlocale("LC_CTYPE")==before)',
  'r <- as.data.frame(readxl::read_xlsx(args[1],sheet="Results",.name_repair="minimal"))',
  'raw <- function(v) vapply(enc2utf8(as.character(v)),function(s) paste(charToRaw(s),collapse=""),character(1),USE.NAMES=FALSE)',
  'stopifnot(identical(raw(names(r)),raw(names(x))), identical(raw(r[[1]]),raw(x[[1]])),',
  '  identical(r[[2]],x[[2]]), identical(r[[3]],x[[3]]))',
  'cat("ok", l10n_info()[["UTF-8"]], "\\n")'), child, useBytes=TRUE)
rscript <- file.path(R.home("bin"),"Rscript")
run <- function(locale, out) {
  # system2(env=...) uses shell-specific assignment on Windows. Inherit a scoped
  # environment instead, and skip user startup files in this isolated test.
  previous_env <- Sys.getenv(c("LC_ALL","LANG"),unset=NA_character_)
  on.exit({
    for(key in names(previous_env)) if(is.na(previous_env[[key]])) Sys.unsetenv(key)
      else do.call(Sys.setenv,setNames(list(previous_env[[key]]),key))
  },add=TRUE)
  Sys.setenv(LC_ALL=locale,LANG=locale)
  log <- suppressWarnings(system2(rscript,c("--vanilla",shQuote(child),shQuote(out),shQuote(locale)),stdout=TRUE,stderr=TRUE,timeout=60))
  status <- attr(log,"status"); if(is.null(status)) status <- 0L
  if(status!=0L) stop("Excel write under ",locale," failed:\n",paste(log,collapse="\n"))
  invisible(log)
}
c_out <- file.path(td,"c.xlsx"); run("C",c_out)
pm_excel_utf8_ctype()
utf8_locale <- Sys.getlocale("LC_CTYPE")
utf8_available <- isTRUE(l10n_info()[["UTF-8"]])
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE","C")))
if(utf8_available) {
  u_out <- file.path(td,"utf8.xlsx"); run(utf8_locale,u_out)
  # Same cells and headers whether R started in C or in UTF-8.
  for(sheet in c("Results","Settings")) {
    a <- readxl::read_xlsx(c_out,sheet=sheet,.name_repair="minimal")
    b <- readxl::read_xlsx(u_out,sheet=sheet,.name_repair="minimal")
    stopifnot(identical(as.data.frame(a),as.data.frame(b)))
  }
}
cat("Excel UTF-8 locale: C-locale headers/cells preserved, LC_CTYPE restored after success and error,",
    "LC_NUMERIC untouched, refusal writes no file, C and UTF-8 workbooks identical.\n")
