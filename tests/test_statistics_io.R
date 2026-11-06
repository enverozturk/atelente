script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
if(.Platform$OS.type=="windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8")))
for(module in c("mapping","report","excel_locale","input_adapter","sample_ratios","normalization_io","statistics","statistics_io"))
  source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
td <- tempfile("pm_statistics_io_");dir.create(td)

peptides <- data.frame(ID=c("001","001E2","10/12","004","005","006"),Sequence=c("PEPTIDE","ACDEFG","KLMNPQ","ACDEF","ACDE","PEPTID"),
  Link="0001",p1=c(0,1,0,3,2,1),p2=c(1,3,4,4,5,2),p3=c(2,4,6,3,7,4),p4=c(4,5,8,5,8,5),check.names=FALSE)
proteins <- data.frame(ID="0001",Members="P1;P1-2",r1=1,r2=1,r3=1,r4=1,check.names=FALSE)
mapping <- data.frame(Sample=c("01","02","03","04"),Condition="",Use="Yes",
  "Peptide column"=paste0("p",1:4),"Protein column"=paste0("r",1:4),check.names=FALSE)
spec <- list(observation_level="peptide_observation",link_separator=NULL,
  peptide=list(id="ID",sequence="Sequence",protein_links="Link",scale="log2",quantity_type="synthetic log2 quantity"),
  protein=list(id="ID",accessions="Members",entity_type="protein_group",scale="linear",quantity_type="synthetic group quantity"))
data <- pm_import_quant_tables(peptides,proteins,spec,mapping);ratios <- pm_sample_ratios(data)
samples <- data.frame(Sample=c("01","02","03","04"),Condition=c("Control","Control","Treatment","Treatment"),Use="Yes",
  "Experimental unit"=c("unit 01","unit 02","unit 03","unit 04"),check.names=FALSE)
comparisons <- data.frame(Name="Treatment vs Control",Treatment="Treatment",Control="Control")
paths <- setNames(file.path(td,c("peptides.csv","samples.csv","comparisons.csv")),c("peptides","samples","comparisons"))
pm_write_csv(peptides,paths[["peptides"]]);pm_write_csv(samples,paths[["samples"]]);pm_write_csv(comparisons,paths[["comparisons"]])
source_hash <- vapply(paths,function(p) digest::digest(file=p,algo="sha256"),character(1))
normalization <- file.path(td,"normalization")
pm_export_sample_ratios(data,ratios,normalization,list(peptides=paths[["peptides"]]),excel=FALSE)
snapshot <- serialize(list(ratios,samples,comparisons),NULL)
result <- pm_test_ratios(ratios$ratios,samples,comparisons,design="independent_units")
stopifnot(identical(snapshot,serialize(list(ratios,samples,comparisons),NULL)),all(result$results$tested))

# Known synthetic export values test precision independently of the fitted p-values.
export_result <- result
export_result$results$p_value[1:2] <- c(1.23456789012345e-120,0.123456789012345)
export_result$results$adjusted_p_value[1:2] <- c(7.4074073407407e-120,0.234567890123456)
export_result$results$prior_df <- rep(NA_real_,nrow(export_result$results))
export_result$results$prior_df_infinite <- rep(TRUE,nrow(export_result$results))
snapshot <- serialize(export_result,NULL)
out <- file.path(td,"excel");csv_only <- file.path(td,"csv_only")
pm_export_statistics(export_result,out,as.list(paths),excel=TRUE,notes="Synthetic export fixture")
pm_export_statistics(export_result,csv_only,as.list(paths),excel=FALSE,notes="Synthetic export fixture")
stopifnot(identical(snapshot,serialize(export_result,NULL)),file.exists(file.path(out,"peptide_statistics.xlsx")),
  !file.exists(file.path(csv_only,"peptide_statistics.xlsx")))
csv_names <- list.files(out,pattern="\\.csv$")
stopifnot(identical(csv_names,sort(c("results.csv","measurements.csv","samples.csv","comparisons.csv"))))
for(name in csv_names) stopifnot(identical(digest::digest(file=file.path(out,name),algo="sha256"),
  digest::digest(file=file.path(csv_only,name),algo="sha256")))
audit <- read.csv(file.path(out,"results.csv"),colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE)
stopifnot(identical(audit$observation_id,export_result$results$observation_id),
  identical(audit$entity_id,export_result$results$entity_id),all(audit$prior_df==""),all(audit$prior_df_infinite=="TRUE"))
near(as.numeric(audit$p_value),export_result$results$p_value)
stopifnot(abs(as.numeric(audit$p_value[1])/export_result$results$p_value[1]-1)<1e-12)
workbook <- file.path(out,"peptide_statistics.xlsx")
sheet <- openxlsx::read.xlsx(workbook,sheet="Results",check.names=FALSE,sep.names=" ")
stopifnot(is.numeric(sheet[["P-value"]]),is.numeric(sheet[["BH-adjusted P-value"]]),
  identical(sheet[["Observation ID"]],export_result$results$observation_id),
  identical(sheet[["Protein / group ID"]],export_result$results$entity_id),
  all(sheet[["Infinite prior degrees of freedom"]]),all(is.na(sheet[["Prior df"]])))
near(sheet[["P-value"]],export_result$results$p_value)
near(sheet[["BH-adjusted P-value"]],export_result$results$adjusted_p_value)
stopifnot(abs(sheet[["P-value"]][1]/export_result$results$p_value[1]-1)<1e-12)
manifest <- jsonlite::read_json(file.path(out,"statistics_manifest.json"),simplifyVector=TRUE)
plain_manifest <- jsonlite::read_json(file.path(csv_only,"statistics_manifest.json"),simplifyVector=TRUE)
stopifnot(manifest$stage=="independent_unit_statistics",manifest$settings$design=="independent_units",
  manifest$notes=="Synthetic export fixture",identical(manifest$notes,plain_manifest$notes),
  manifest$settings$limma_version==as.character(utils::packageVersion("limma")),
  manifest$row_counts$results==6L,manifest$row_counts$measurements==24L,
  manifest$workbook_columns$Results$p_value=="P-value",manifest$excel_export,!plain_manifest$excel_export)
for(name in names(paths)) stopifnot(manifest$inputs[[name]]$sha256==source_hash[[name]],
  manifest$inputs[[name]]$file==normalizePath(paths[[name]],winslash="/"))
for(name in names(manifest$csv_sha256)) stopifnot(manifest$csv_sha256[[name]]==
  digest::digest(file=file.path(out,paste0(name,".csv")),algo="sha256"))
protected <- digest::digest(file=file.path(out,"results.csv"),algo="sha256")
fails(pm_export_statistics(result,out,as.list(paths),excel=FALSE))
stopifnot(protected==digest::digest(file=file.path(out,"results.csv"),algo="sha256"))
fails(pm_export_statistics(result,file.path(td,"invalid_excel"),as.list(paths),excel="false"))
bad <- export_result;bad$results$prior_df[1] <- Inf
fails(pm_export_statistics(bad,file.path(td,"unencoded_inf"),as.list(paths),excel=FALSE))
stopifnot(!dir.exists(file.path(td,"unencoded_inf")))

# Hash verification precedes parsing, and signed typed values/literal IDs survive reload.
loaded <- pm_load_ratio_statistics(normalization)
stopifnot(identical(loaded$ratios$observation_id,ratios$ratios$observation_id),
  identical(loaded$ratios$entity_id,ratios$ratios$entity_id),identical(loaded$ratios$sample_id,ratios$ratios$sample_id),
  is.numeric(loaded$ratios$log2_ratio),is.logical(loaded$ratios$usable_pair),
  all(loaded$ratios$condition==""),!anyNA(loaded$ratios$condition))
near(loaded$ratios$log2_ratio,ratios$ratios$log2_ratio)
near(pm_test_ratios(loaded$ratios,samples,comparisons)$results$p_value,result$results$p_value)
tamper <- file.path(td,"tamper");dir.create(tamper)
stopifnot(all(file.copy(file.path(normalization,c("sample_ratios.csv","normalization_manifest.json")),tamper)))
cat("\n",file=file.path(tamper,"sample_ratios.csv"),append=TRUE)
fails(pm_load_ratio_statistics(tamper))
bad_manifest <- jsonlite::read_json(file.path(normalization,"normalization_manifest.json"),simplifyVector=TRUE)
bad_ratios <- ratios$ratios;bad_ratios$log2_ratio <- as.character(bad_ratios$log2_ratio);bad_ratios$log2_ratio[1] <- "not numeric"
pm_write_csv(bad_ratios,file.path(tamper,"sample_ratios.csv"))
bad_manifest$csv_sha256$sample_ratios <- digest::digest(file=file.path(tamper,"sample_ratios.csv"),algo="sha256")
jsonlite::write_json(bad_manifest,file.path(tamper,"normalization_manifest.json"),auto_unbox=TRUE,null="null")
fails(pm_load_ratio_statistics(tamper))
bad_ratios <- ratios$ratios;bad_ratios$usable_pair <- "yes";pm_write_csv(bad_ratios,file.path(tamper,"sample_ratios.csv"))
bad_manifest$csv_sha256$sample_ratios <- digest::digest(file=file.path(tamper,"sample_ratios.csv"),algo="sha256")
jsonlite::write_json(bad_manifest,file.path(tamper,"normalization_manifest.json"),auto_unbox=TRUE,null="null")
fails(pm_load_ratio_statistics(tamper));fails(pm_load_ratio_statistics(file.path(td,"missing")))

# CLI requires an explicit contrast and design, using configuration-relative paths.
config_path <- file.path(td,"statistics.json")
config <- list(ratios_dir="normalization",samples="samples.csv",comparisons="comparisons.csv",output_dir="cli_results",
  design="independent_units",excel=FALSE,notes="Synthetic CLI test")
save_config <- function(x) jsonlite::write_json(x,config_path,auto_unbox=TRUE,pretty=TRUE,null="null")
run_cli <- function() {
  executable <- if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript"
  output <- suppressWarnings(system2(file.path(R.home("bin"),executable),
    c(shQuote(file.path(root,"statistics.R")),"--config",shQuote(config_path)),stdout=TRUE,stderr=TRUE))
  status <- attr(output,"status");if(is.null(status)) status <- 0L
  list(status=status,output=output)
}
save_config(config);cli <- run_cli()
if(cli$status!=0L) stop(paste(cli$output,collapse="\n"))
cli_results <- read.csv(file.path(td,"cli_results","results.csv"),colClasses="character",fileEncoding="UTF-8-BOM",check.names=FALSE)
near(as.numeric(cli_results$p_value),result$results$p_value)
stopifnot(identical(cli_results$observation_id,result$results$observation_id),all(cli_results$status=="tested"))
cm <- jsonlite::read_json(file.path(td,"cli_results","statistics_manifest.json"),simplifyVector=TRUE)
stopifnot(cm$inputs$configuration$sha256==digest::digest(file=config_path,algo="sha256"),
  cm$inputs$ratio_manifest$sha256==digest::digest(file=file.path(normalization,"normalization_manifest.json"),algo="sha256"),
  !cm$excel_export,run_cli()$status!=0L)
for(design in c("paired","batch")) {
  bad <- config;bad$output_dir <- paste0("invalid_",design);bad$design <- design;save_config(bad)
  stopifnot(run_cli()$status!=0L,!dir.exists(file.path(td,bad$output_dir)))
}
bad <- config;bad$output_dir <- "invalid_excel_cli";bad$excel <- "false";save_config(bad)
stopifnot(run_cli()$status!=0L,!dir.exists(file.path(td,bad$output_dir)))
bad <- config;bad$comparisons <- NULL;save_config(bad);stopifnot(run_cli()$status!=0L)
bad <- config;bad$design <- NULL;save_config(bad);stopifnot(run_cli()$status!=0L)
stopifnot(identical(source_hash,vapply(paths,function(p) digest::digest(file=p,algo="sha256"),character(1))))
stopifnot(startsWith(normalizePath(td,winslash="/"),paste0(normalizePath(tempdir(),winslash="/"),"/")))
unlink(td,recursive=TRUE)
cat("Statistics IO tests passed.\n")
