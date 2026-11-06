script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
source(file.path(root,"R/mapping.R"))
td <- tempfile("csv_delimiters_");dir.create(td)
x <- data.frame(Sequence=c("ACDEK","MPEPTIDER"), Protein=c("P00001;P00002","P00003"),
  Description=c('quoted, comma; and "quote"',"two\nlines"), Quantity=c("100.5","0"),check.names=FALSE)
for(sep in c(",",";")) {
  path<-file.path(td,if(sep==",") "comma.csv" else "semicolon.csv")
  write.table(x,path,sep=sep,row.names=FALSE,quote=TRUE,qmethod="double",fileEncoding="UTF-8")
  y<-pm_read_table(path)
  stopifnot(identical(names(x),names(y)),identical(as.matrix(x),as.matrix(y)))
}
# An unquoted semicolon header and decimal-comma cell are distinct conventions.
path<-file.path(td,"decimal.csv")
writeLines(c("Sequence;Protein;Quantity","ACDEK;P00001;100,5"),path)
y<-pm_read_table(path)
stopifnot(identical(names(y),c("Sequence","Protein","Quantity")),y$Quantity=="100,5")
err<-tryCatch(pm_validate_quantity_columns(y,"Quantity"),error=identity)
stopifnot(inherits(err,"error"),grepl("decimal point",conditionMessage(err)))
# A delimiter in a quoted header is literal, and a single-column CSV remains valid.
writeLines(c('"Description; detail",Quantity','"a;b",1.5'),path)
stopifnot(identical(names(pm_read_table(path)),c("Description; detail","Quantity")))
writeLines(c("Sequence","ACDEK"),path);stopifnot(pm_read_table(path)$Sequence=="ACDEK")
writeLines(c('"Sequence', 'header";Quantity','ACDEK;1.5'),path)
stopifnot(identical(names(pm_read_table(path)),c("Sequence\nheader","Quantity")))
# A sampled data row ending inside a quoted multiline value must not affect
# header detection or emit a spurious EOF-in-quote warning.
writeLines(c("Sequence;Description",rep("ACDEK;plain",61L),'ACDEK;"first','second','third"'),path)
old_warn<-getOption("warn");options(warn=2)
stopifnot(nrow(pm_read_table(path))==62L);options(warn=old_warn)
writeLines(c("Sequence,Protein;Quantity","ACDEK,P00001;1.5"),path)
err<-tryCatch(pm_read_table(path),error=identity)
stopifnot(inherits(err,"error"),grepl("Ambiguous CSV delimiter",conditionMessage(err)))
cat("Quoted comma/semicolon CSV, literal fields, decimal validation and ambiguity checks passed.\n")
