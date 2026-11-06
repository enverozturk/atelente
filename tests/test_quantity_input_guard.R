# Full workflow rejects malformed selected quantities before producing a result.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash="/")))
api <- new.env(parent=globalenv())
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))sys.source(f,api)
td <- tempfile("quantity_guard_")
cfg <- api$pm_wizard_demo(td)
src <- api$pm_read_table(file.path(td,"peptides.csv"))
check_error <- function(values, scale="linear", native=FALSE) {
  input <- src;input[["Signal 1"]] <- values
  run <- cfg;run$column_mapping$quantity_scale <- scale
  column <- "Signal 1"
  if(native) {
    input <- data.frame("Annotated Sequence"=src[["Peptide sequence"]],
      "Master Protein Accessions"=src[["Representative ID"]],
      "Abundance: sample"=values,check.names=FALSE)
    run$column_mapping <- NULL;run$sample_design <- NULL;run$assignment_mode <- "unique"
    column <- "Abundance: sample"
  }
  api$pm_write_csv(input,file.path(td,"peptides.csv"))
  before <- digest::digest(file=file.path(td,"peptides.csv"),algo="sha256")
  err <- tryCatch({api$pm_run(run,td,file.path(root,"inst","templates","report_template.html"));NULL},error=identity)
  stopifnot(inherits(err,"error"))
  msg <- conditionMessage(err)
  stopifnot(grepl(column,msg,fixed=TRUE),grepl("row",msg,fixed=TRUE),
    grepl("decimal point",msg,fixed=TRUE),!dir.exists(file.path(td,"results")),
    identical(before,digest::digest(file=file.path(td,"peptides.csv"),algo="sha256")))
}
check_error(rep("100,5",6))
check_error(c("100.5","100,5","20","30",NA,"0"))
check_error(c("1","2,5","3","4","5","6"),scale="log2")
check_error(c("100.5","100,5","20","30",NA,"0"),native=TRUE)
check_error(c("100.5","oops","20","30",NA,"0"))
# Excluded columns may retain diagnostic input, but selected raw values remain exact.
input <- src;input[["Signal 1"]] <- c("0",NA,"1.25e-12","23.456789012345","-2","100")
input[["Unused"]] <- "1,25"
api$pm_write_csv(input,file.path(td,"peptides.csv"))
api$pm_run(cfg,td,file.path(root,"inst","templates","report_template.html"))
m <- api$pm_read_table(file.path(td,"results","audit","measurements.csv"))
actual <- as.numeric(m$value[m$source_column=="Abundance: Selected_001"])
stopifnot(isTRUE(all.equal(actual,c(0,NA,1.25e-12,23.456789012345,-2,100),tolerance=1e-14)))
cat("Selected malformed quantities fail clearly; valid/missing/zero quantities remain unchanged.\n")
