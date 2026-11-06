# UTF-8 character handling for Excel writing.
#
# openxlsx escapes non-ASCII text in headers and cells (as <U+....> or <xx>
# byte codes) when the session LC_CTYPE is not UTF-8, e.g. LC_ALL=C on Linux,
# Docker or CI. Text is converted while cells are written, so the switch must
# cover workbook creation, writeData() and saveWorkbook(), not only the save.
#
# Contract:
# - An existing UTF-8 session is left untouched.
# - Only LC_CTYPE changes; LC_NUMERIC and other categories are not modified.
# - A candidate counts only if UTF-8 is actually active after switching.
# - The previous LC_CTYPE is restored on success and on error.
# - Without a usable UTF-8 locale, stop before any output file is created.
# This is a per-call session setting, not a persistent application start-up option.

pm_excel_utf8_candidates <- c("C.UTF-8","C.utf8","en_US.UTF-8","en_US.utf8",".UTF-8")

# Switch LC_CTYPE to UTF-8 if needed. Returns the previous LC_CTYPE when it was
# changed, or NULL when the session was already UTF-8. Callers must restore it.
pm_excel_utf8_ctype <- function(candidates=pm_excel_utf8_candidates) {
  if(isTRUE(l10n_info()[["UTF-8"]])) return(NULL)
  old <- Sys.getlocale("LC_CTYPE")
  for(candidate in candidates)
    if(nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE",candidate))) && isTRUE(l10n_info()[["UTF-8"]]))
      return(old)
  suppressWarnings(Sys.setlocale("LC_CTYPE",old))
  stop("Excel export requires a UTF-8 locale to preserve non-ASCII text, but none is available. ",
       "Start R with a UTF-8 locale (for example LANG=en_US.UTF-8) or disable Excel export.",call.=FALSE)
}

# Evaluate Excel-writing code (createWorkbook ... writeData ... saveWorkbook)
# with a UTF-8 LC_CTYPE, restoring the caller's LC_CTYPE afterwards.
# `code` is evaluated lazily, i.e. only after the locale switch succeeded, so an
# unavailable UTF-8 locale fails before any file is opened or written.
pm_with_excel_utf8 <- function(code,candidates=pm_excel_utf8_candidates) {
  previous <- pm_excel_utf8_ctype(candidates)
  if(!is.null(previous)) on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE",previous)),add=TRUE)
  code
}
