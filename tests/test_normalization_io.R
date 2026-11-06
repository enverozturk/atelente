script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
if(.Platform$OS.type=="windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8")))
for(module in c("mapping","report","excel_locale","input_adapter","sample_ratios","normalization_io"))
  source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
td <- tempfile("pm_normalization_io_");dir.create(td)

peptides <- data.frame(ID=c("001","001E2","10/12","missing-link"),Sequence=c("PEPTIDE","PEPTIDE","ACDEFG","KLMNPQ"),
  Links=c("0001","0001","0002","unknown"),Mods=c("","Oxidation [M1]","",""),
  "log2 replicate A"=c("0","-2","3","1"),"log2 replicate B"=c("4","","oops","-1"),
  Note=c("\u00e7\u0131\u011f","comma, quote\"","keep","all"),check.names=FALSE)
proteins <- data.frame(Group=c("0001","0002"),Members=c("P1;P1-2","P2"),
  "intensity A"=c("2","0"),"intensity B"=c("4","8"),check.names=FALSE)
samples <- data.frame(Sample=c("01","02"),Condition="same condition",Use="Yes",
  "Peptide column"=c("log2 replicate A","log2 replicate B"),
  "Protein column"=c("intensity A","intensity B"),check.names=FALSE)
spec <- list(observation_level="peptide_observation",link_separator=NULL,
  peptide=list(id="ID",sequence="Sequence",protein_links="Links",modifications="Mods",scale="log2",quantity_type="reported peptide quantity"),
  protein=list(id="Group",accessions="Members",entity_type="protein_group",scale="linear",quantity_type="reported group intensity"))
paths <- setNames(file.path(td,c("peptides.csv","proteins.csv","samples.csv")),c("peptides","proteins","samples"))
pm_write_csv(peptides,paths[["peptides"]]);pm_write_csv(proteins,paths[["proteins"]]);pm_write_csv(samples,paths[["samples"]])
data <- pm_import_quant_tables(peptides,proteins,spec,samples);result <- pm_sample_ratios(data)
before <- serialize(list(data,result),NULL)
out <- file.path(td,"with_excel");csv_only <- file.path(td,"csv_only")
manifest <- pm_export_sample_ratios(data,result,out,as.list(paths),excel=TRUE,notes="Synthetic IO test")
pm_export_sample_ratios(data,result,csv_only,as.list(paths),excel=FALSE,notes="Synthetic IO test")
stopifnot(identical(before,serialize(list(data,result),NULL)),file.exists(file.path(out,"sample_ratios.xlsx")),
  !file.exists(file.path(csv_only,"sample_ratios.xlsx")))
csv_names <- list.files(out,pattern="\\.csv$")
stopifnot(length(csv_names)==10L,identical(csv_names,list.files(csv_only,pattern="\\.csv$")))
for(name in csv_names) stopifnot(identical(digest::digest(file=file.path(out,name),algo="sha256"),
  digest::digest(file=file.path(csv_only,name),algo="sha256")))
audit <- pm_read_table(file.path(out,"sample_ratios.csv"))
stopifnot(identical(audit$observation_id,result$ratios$observation_id),
  identical(audit$entity_id,result$ratios$entity_id),identical(audit$sample_id,result$ratios$sample_id),
  identical(audit$ratio_status,result$ratios$ratio_status),nrow(audit)==8L)
near(as.numeric(audit$normalized_ratio),result$ratios$normalized_ratio)
source_csv <- pm_read_table(file.path(out,"source_peptides.csv"))
stopifnot(identical(source_csv$ID,peptides$ID),identical(source_csv$Note,peptides$Note),
  source_csv[["log2 replicate A"]][1]=="0",is.na(source_csv[["log2 replicate B"]][2]),
  source_csv[["log2 replicate B"]][3]=="oops")
measurements <- pm_read_table(file.path(out,"peptide_measurements.csv"))
stopifnot(measurements$original_value[1]=="0",measurements$value_status[1]=="observed",
  measurements$value_status[6]=="missing",measurements$value_status[7]=="invalid_numeric_value")

# Excel cells remain numeric for quantities and literal text for identifiers.
workbook <- file.path(out,"sample_ratios.xlsx")
stopifnot(identical(openxlsx::getSheetNames(workbook),
  c("Sample ratios","Peptide measurements","Protein measurements","Samples","Read me")))
sheet <- openxlsx::read.xlsx(workbook,sheet="Sample ratios",check.names=FALSE,sep.names=" ")
stopifnot(is.numeric(sheet[["Peptide / protein"]]),is.numeric(sheet[["Log2 ratio"]]),
  identical(sheet[["Observation ID"]],result$ratios$observation_id),
  identical(sheet[["Protein / group ID"]],result$ratios$entity_id),
  identical(sheet$Sample,result$ratios$sample_id))
near(sheet[["Peptide / protein"]],result$ratios$normalized_ratio)
near(sheet[["Log2 ratio"]],result$ratios$log2_ratio)
near(sheet[["Peptide / protein"]][1:4],c(0.5,4,0.125,NA))
pm <- openxlsx::read.xlsx(workbook,sheet="Peptide measurements",check.names=FALSE,sep.names=" ")
stopifnot(is.numeric(pm[["Input quantity"]]),is.numeric(pm[["Log2 quantity"]]),
  pm[["Original value text"]][1]=="0",is.na(pm[["Input quantity"]][6]),
  pm[["Original value text"]][7]=="oops")
near(pm[["Log2 quantity"]][1:2],c(0,-2))

# Manifest records exact scale/settings, literal workbook labels, row counts and hashes.
m <- jsonlite::read_json(file.path(out,"normalization_manifest.json"),simplifyVector=TRUE)
stopifnot(m$stage=="per_sample_ratios",m$spec$peptide$scale=="log2",m$spec$protein$scale=="linear",
  "link_separator" %in% names(m$spec),is.null(m$spec$link_separator),
  m$settings$scales$peptide=="log2",m$settings$scales$protein=="linear",
  m$settings$operation==result$settings$operation,m$excel_export,
  m$row_counts$sample_ratios==nrow(result$ratios),m$row_counts$entities==2L,
  m$workbook_columns[["Sample ratios"]]$normalized_ratio=="Peptide / protein",
  m$notes=="Synthetic IO test")
for(name in names(paths)) stopifnot(m$inputs[[name]]$file==normalizePath(paths[[name]],winslash="/"),
  m$inputs[[name]]$sha256==digest::digest(file=paths[[name]],algo="sha256"))
for(name in names(m$csv_sha256)) stopifnot(m$csv_sha256[[name]]==digest::digest(file=file.path(out,paste0(name,".csv")),algo="sha256"))
stopifnot(!jsonlite::read_json(file.path(csv_only,"normalization_manifest.json"))$excel_export)
protected <- digest::digest(file=file.path(out,"sample_ratios.csv"),algo="sha256")
fails(pm_export_sample_ratios(data,result,out,as.list(paths),excel=FALSE))
stopifnot(protected==digest::digest(file=file.path(out,"sample_ratios.csv"),algo="sha256"))
fails(pm_export_sample_ratios(data,result,file.path(td,"invalid_excel"),as.list(paths),excel="false"))
stopifnot(!dir.exists(file.path(td,"invalid_excel")))

# The CLI resolves configuration-relative files and accepts no contrast or FASTA requirement.
config_path <- file.path(td,"normalization.json")
config <- list(peptides="peptides.csv",proteins="proteins.csv",samples="samples.csv",output_dir="cli_results",
  spec=spec,excel=FALSE,notes="Synthetic CLI test")
save_config <- function(x) jsonlite::write_json(x,config_path,auto_unbox=TRUE,pretty=TRUE,null="null")
run_cli <- function() {
  executable <- if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript"
  output <- suppressWarnings(system2(file.path(R.home("bin"),executable),
    c(shQuote(file.path(root,"normalize.R")),"--config",shQuote(config_path)),stdout=TRUE,stderr=TRUE))
  status <- attr(output,"status");if(is.null(status)) status <- 0L
  list(status=status,output=output)
}
save_config(config);cli <- run_cli()
# Diagnostic low-level imports above remain available, but the full CLI cannot
# silently discard malformed cells and publish a successful analysis.
stopifnot(cli$status!=0L,any(grepl("decimal point",cli$output,fixed=TRUE)),
  !dir.exists(file.path(td,"cli_results")))
peptides[["log2 replicate B"]][3] <- NA_character_
pm_write_csv(peptides,paths[["peptides"]])
cli <- run_cli()
if(cli$status!=0L) stop(paste(cli$output,collapse="\n"))
cm <- jsonlite::read_json(file.path(td,"cli_results","normalization_manifest.json"),simplifyVector=TRUE)
stopifnot(cm$inputs$configuration$sha256==digest::digest(file=config_path,algo="sha256"),
  cm$row_counts$sample_ratios==8L,!cm$excel_export)
stopifnot(run_cli()$status!=0L)
bad <- config;bad$output_dir <- "bad_config_output";bad$excel <- "false";save_config(bad)
stopifnot(run_cli()$status!=0L,!dir.exists(file.path(td,"bad_config_output")))
bad <- config;bad$output_dir <- "bad_mapping_output";bad$spec$peptide$sequence <- "absent";save_config(bad)
stopifnot(run_cli()$status!=0L,!dir.exists(file.path(td,"bad_mapping_output")))
bad <- config;bad$proteins <- NULL;save_config(bad);stopifnot(run_cli()$status!=0L)

stopifnot(startsWith(normalizePath(td,winslash="/"),paste0(normalizePath(tempdir(),winslash="/"),"/")))
unlink(td,recursive=TRUE)
cat("Normalization IO tests passed.\n")
