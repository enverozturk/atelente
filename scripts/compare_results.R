# Compare two Atelente result workbooks from the same input and settings.
#
#   Rscript scripts/compare_results.R old/atelente_results.xlsx new/atelente_results.xlsx
#
# Checks that every sheet except Settings has the same rows and that every
# column of the old workbook is unchanged in the new one. Columns that exist only
# in the new workbook are listed. When the new workbook has
# "Spanning terminal form", its counts are summarised for the Results sheet.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript scripts/compare_results.R OLD.xlsx NEW.xlsx")
old_path <- normalizePath(args[1],mustWork=TRUE); new_path <- normalizePath(args[2],mustWork=TRUE)
read_sheet <- function(path,sheet) openxlsx::read.xlsx(path,sheet=sheet,check.names=FALSE,sep.names=" ",na.strings=character())
old_sheets <- openxlsx::getSheetNames(old_path); new_sheets <- openxlsx::getSheetNames(new_path)
problems <- character()
if(!identical(old_sheets,new_sheets)) problems <- c(problems,sprintf("Sheets differ: old [%s], new [%s]",
  paste(old_sheets,collapse=", "),paste(new_sheets,collapse=", ")))
for(sheet in setdiff(intersect(old_sheets,new_sheets),"Settings")){
  o <- read_sheet(old_path,sheet); n <- read_sheet(new_path,sheet)
  added <- setdiff(names(n),names(o)); removed <- setdiff(names(o),names(n))
  cat(sprintf("\n[%s] rows old %d / new %d\n",sheet,nrow(o),nrow(n)))
  if(length(added)) cat("  Added columns:",paste(added,collapse=", "),"\n")
  if(length(removed)) problems <- c(problems,sprintf("[%s] removed columns: %s",sheet,paste(removed,collapse=", ")))
  if(nrow(o)!=nrow(n)) {problems <- c(problems,sprintf("[%s] row count differs",sheet));next}
  for(col in intersect(names(o),names(n))){
    a <- o[[col]]; b <- n[[col]]
    same <- if(is.numeric(a)&&is.numeric(b)) isTRUE(all.equal(a,b,tolerance=1e-12)) else identical(as.character(a),as.character(b))
    if(!same) problems <- c(problems,sprintf("[%s] column changed: %s (%d rows differ)",sheet,col,
      sum(xor(is.na(a),is.na(b))|(!is.na(a)&!is.na(b)&as.character(a)!=as.character(b)))))
  }
  if(!length(grep(sprintf("^\\[%s\\]",sheet),problems))) cat("  All old columns unchanged.\n")
}
if("Results" %in% new_sheets){
  r <- read_sheet(new_path,"Results")
  if("Spanning terminal form" %in% names(r)){
    cat("\nSpanning terminal form in Results (all rows):\n"); print(table(Role=r$Role,`Spanning terminal form`=r$`Spanning terminal form`,useNA="ifany"))
    if(all(c("Test status","BH-adjusted P") %in% names(r))){
      tested <- r$`Test status` %in% "Tested"
      cat("\nTested representative pairs:\n"); print(table(Role=r$Role[tested],`Spanning terminal form`=r$`Spanning terminal form`[tested],useNA="ifany"))
      sig <- tested & r$`BH-adjusted P` < 0.05
      cat("\nTested pairs with q < 0.05:\n"); print(table(Role=r$Role[sig],`Spanning terminal form`=r$`Spanning terminal form`[sig],useNA="ifany"))
    }
  }
}
if(length(problems)){
  cat("\nDIFFERENCES FOUND:\n",paste0("  ",problems,"\n"),sep=""); quit(status=1)
}
cat("\nOK: results are unchanged apart from added columns.\n")
