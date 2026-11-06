script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
source(file.path(root,"R","report.R"),encoding="UTF-8")
td <- tempfile("pm_maxquant_cli_");dir.create(td)
peptides <- data.frame(Sequence=c("PEPTIDE","PEPTIDER"),id=c("0","01"),
  "Protein group IDs"=c("001","002"),"Potential contaminant"=c("","+"),
  "Intensity A"=c("100","200"),"Intensity B"=c("200","0"),check.names=FALSE)
proteins <- data.frame(id=c("001","002"),"Protein IDs"=c("P1;P2","P3"),
  "Intensity A"=c("1000","1000"),"Intensity B"=c("2000","2000"),check.names=FALSE)
samples <- data.frame(Sample=c("A","B"),Condition="",Use="Yes",
  "Peptide column"=c("Intensity A","Intensity B"),"Protein column"=c("Intensity A","Intensity B"),check.names=FALSE)
for(name in c("peptides","proteins","samples")) pm_write_csv(get(name),file.path(td,paste0(name,".csv")))
config <- list(adapter="maxquant_peptides",peptides="peptides.csv",proteins="proteins.csv",samples="samples.csv",
  output_dir="results",excel=FALSE,maxquant_options=list(protein_quantity="Intensity",exclude_only_identified_by_site=TRUE))
run <- function(cfg) {
  jsonlite::write_json(cfg,file.path(td,"config.json"),auto_unbox=TRUE,pretty=TRUE)
  status <- suppressWarnings(system2(file.path(R.home("bin"),"Rscript"),
    c(shQuote(file.path(root,"normalize.R")),"--config",shQuote(file.path(td,"config.json"))),
    stdout=file.path(td,"stdout.txt"),stderr=file.path(td,"stderr.txt")))
  status
}
status <- run(config)
if(status!=0L) stop(paste(readLines(file.path(td,"stderr.txt")),collapse="\n"))
out <- read.csv(file.path(td,"results","sample_ratios.csv"),colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM")
stopifnot(identical(out$observation_id,c("0","0","01","01")),identical(out$entity_id,c("001","001","002","002")),
  all(abs(as.numeric(out$normalized_ratio[1:2])-.1)<1e-12),all(out$ratio_status[3:4]=="excluded_observation"))
original <- read.csv(file.path(td,"results","source_peptides.csv"),colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM")
stopifnot(identical(original[["Potential contaminant"]],c("","+")))
manifest <- jsonlite::read_json(file.path(td,"results","normalization_manifest.json"),simplifyVector=TRUE)
stopifnot(manifest$spec$protein$entity_type=="protein_group",manifest$spec$peptide$scale=="linear")
bad <- config;bad$output_dir <- "invalid";bad$adapter <- "guess";stopifnot(run(bad)!=0L)
bad <- config;bad$output_dir <- "invalid";bad$spec <- list(dummy=TRUE);stopifnot(run(bad)!=0L)
bad <- config;bad$output_dir <- "invalid";bad$maxquant_options <- list(guess=TRUE);stopifnot(run(bad)!=0L)
cat("MaxQuant CLI tests passed.\n")
