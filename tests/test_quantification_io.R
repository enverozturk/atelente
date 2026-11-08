script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
for(module in c("mapping","report","excel_locale","excel","sample_metadata","quantification","quantification_io")) source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
td <- tempfile("pm_quant_io_");dir.create(td)
fasta <- file.path(td,"test.fasta");writeLines(c(">sp|P1|TEST","MAAAAPEPTIDER"),fasta)
peptides <- data.frame("Annotated Sequence"=c("AAA","PEPTIDE"),"Master Protein Accessions"="P1",
  "Abundance: F1: Sample"=c("10","20"),"Abundance: F2: Sample"=c("40","0"),check.names=FALSE)
proteins <- data.frame(Accession="P1","Abundance: F1: Sample"="100","Abundance: F2: Sample"="200",check.names=FALSE)
mapped <- pm_map(peptides,pm_read_fasta(fasta),proteins)
export <- file.path(td,"mapping")
pm_export(mapped,export,file.path(root,"inst","templates","report_template.html"),list(fasta=fasta),excel=TRUE)
stopifnot(file.exists(file.path(export,"sample_metadata.xlsx")))
loaded <- pm_load_quantification_result(export)
stopifnot(identical(mapped$measurements$value,loaded$measurements$value),identical(mapped$matches$is_master,loaded$matches$is_master),
  identical(as.character(mapped$source_proteins$Accession),loaded$source_proteins$Accession))
metadata <- pm_sample_template(loaded)
metadata$samples$Condition <- c("Control","Treatment")
metadata$comparisons <- data.frame(Name="TvsC",Treatment="Treatment",Control="Control")
sample_path <- file.path(td,"samples.xlsx");pm_write_sample_metadata(metadata,sample_path)
metadata <- pm_read_sample_metadata(sample_path,loaded,TRUE)
q <- pm_adjust_peptide_effects(loaded,metadata)
stopifnot(isTRUE(all.equal(q$effects$adjusted_log2fc[q$effects$sequence=="AAA"],1,tolerance=1e-12)),
  is.na(q$effects$adjusted_log2fc[q$effects$sequence=="PEPTIDE"]))
out <- file.path(td,"quantification")
pm_export_quantification(q,metadata,out,list(samples=sample_path),assumptions="Synthetic test",protein_abundance_source="Synthetic supplied protein table")
sheet <- openxlsx::read.xlsx(file.path(out,"protein_adjusted_effects.xlsx"),sheet="Peptide effects",check.names=FALSE,sep.names=" ")
stopifnot(isTRUE(all.equal(sheet[["Adjusted log2FC"]],q$effects$adjusted_log2fc,tolerance=1e-12)),
  file.exists(file.path(out,"quantification_manifest.json")))
# Missing protein columns and inconsistent numeric exports cannot silently run.
bad <- file.path(td,"bad");dir.create(bad)
file.copy(list.files(export,pattern="\\.csv$",full.names=TRUE),bad)
m <- loaded$measurements;m$value <- as.character(m$value);m$value[1] <- "bad-value"
pm_write_csv(m,file.path(bad,"measurements.csv"))
stopifnot(inherits(try(pm_load_quantification_result(bad),silent=TRUE),"try-error"))
# Run the real CLI against flat and compact layouts. Its manifest must hash the
# files actually consumed; no mixture or silent fallback to another export.
executable <- file.path(R.home("bin"),"Rscript")
invoke <- function(directory,destination) {
  config <- file.path(td,"quantification.json")
  jsonlite::write_json(list(results_dir=directory,samples=sample_path,output_dir=destination),config,auto_unbox=TRUE)
  log <- suppressWarnings(system2(executable,shQuote(c(file.path(root,"scripts","quantify.R"),"--config",config)),stdout=TRUE,stderr=TRUE))
  status <- attr(log,"status");if(is.null(status))status<-0L
  list(status=status,log=paste(log,collapse="\n"))
}
sha <- function(path) digest::digest(file=path,algo="sha256")
flat_out <- file.path(td,"flat_cli");execution<-invoke(export,flat_out)
if(execution$status!=0L)stop(execution$log)
audit <- file.path(export,"audit");dir.create(audit)
csvs <- list.files(export,pattern="\\.csv$",full.names=TRUE)
stopifnot(all(file.rename(csvs,file.path(audit,basename(csvs)))))
compact_loaded <- pm_load_quantification_result(export)
without_path <- function(result) lapply(result,function(table){attr(table,"source_file")<-NULL;table})
stopifnot(identical(without_path(compact_loaded),without_path(loaded)))
compact_out <- file.path(td,"compact_cli");execution<-invoke(export,compact_out)
if(execution$status!=0L)stop(execution$log)
for(file in c("effects.csv","measurement_pairs.csv","samples.csv","comparisons.csv"))
  stopifnot(sha(file.path(flat_out,file))==sha(file.path(compact_out,file)))
for(sheet in openxlsx::getSheetNames(file.path(flat_out,"protein_adjusted_effects.xlsx")))
  stopifnot(identical(readxl::read_excel(file.path(flat_out,"protein_adjusted_effects.xlsx"),sheet=sheet),
    readxl::read_excel(file.path(compact_out,"protein_adjusted_effects.xlsx"),sheet=sheet)))
manifest <- jsonlite::read_json(file.path(compact_out,"quantification_manifest.json"))
for(name in c("peptides","matches","measurements","proteins")) {
  input <- manifest$inputs[[name]]
  stopifnot(dirname(input$file)==normalizePath(audit,winslash="/"),sha(input$file)==input$sha256)
}
stopifnot(identical(without_path(pm_load_quantification_result(audit)),without_path(loaded)))
# A second layout is ambiguous even when its file is an identical copy.
stopifnot(file.copy(file.path(audit,"peptides.csv"),file.path(export,"peptides.csv")))
refused <- file.path(td,"ambiguous_cli");execution<-invoke(export,refused)
stopifnot(execution$status!=0L,grepl("Ambiguous mapping exports",execution$log,fixed=TRUE),!dir.exists(refused))
unlink(file.path(export,"peptides.csv"))
stopifnot(file.rename(file.path(audit,"source_proteins.csv"),file.path(audit,"saved_proteins.csv")))
refused <- file.path(td,"missing_protein_cli");execution<-invoke(export,refused)
stopifnot(execution$status!=0L,grepl("source_proteins.csv",execution$log,fixed=TRUE),!dir.exists(refused))
stopifnot(startsWith(normalizePath(td,winslash="/"),paste0(normalizePath(tempdir(),winslash="/"),"/")))
unlink(td,recursive=TRUE)
cat("Quantification IO tests passed.\n")
