script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))source(f,encoding="UTF-8")
template <- file.path(root,"inst","templates","report_template.html")
# A spanning form may itself be a terminal form. The Results sheet reports its
# ends without changing selection, ratios or tests.
run_case <- function(fasta,peptides,cleavage){
  td <- tempfile("spanning_ends_");pm_wizard_demo(td)
  cfg <- jsonlite::read_json(file.path(td,"mapping.json"),simplifyVector=TRUE)
  writeLines(fasta,file.path(td,"database.fasta"))
  pm_write_csv(peptides,file.path(td,"peptides.csv"))
  cfg$excel <- TRUE;cfg$cleavage <- cleavage
  pm_run(cfg,config_base=td,template=template)
  r <- openxlsx::read.xlsx(file.path(td,"results","atelente_results.xlsx"),sheet="Results",check.names=FALSE,sep.names=" ")
  r$key <- paste(r$Bond,r$Role,r$`Neo sequence`)
  r
}
quant <- function(control,treatment){
  q <- cbind(outer(control,c(1,1.08,0.93)),outer(treatment,c(1.05,0.95,1.02)))
  colnames(q) <- paste0("Signal ",1:6);as.data.frame(q)
}
row_for <- function(r,key){i <- match(key,r$key);stopifnot(!is.na(i));r[i,]}

# Unlabelled tryptic data. TEST1 has an intact tryptic spanning form; TEST2 has
# only the longer neo-N form, which therefore remains the available spanning form.
seqs <- c("FGAPSASVAAAPAPAR","GAPSASVAAAPAPAR","APSASVAAAPAPAR","GLPDTNVQEHPYQGR","LPDTNVQEHPYQGR")
ids <- c("TEST1","TEST1","TEST1","TEST2","TEST2")
p <- data.frame("Peptide sequence"=seqs,"Modification detail"="","Representative ID"=ids,"Reported IDs"=ids,check.names=FALSE)
p <- cbind(p,quant(c(1000,100,100,100,100),c(250,400,400,400,400)))
r <- run_case(c(">sp|TEST1|ONE GN=ONE","MAGKFGAPSASVAAAPAPARLLSEK",">sp|TEST2|TWO GN=TWO","MSTKWGLPDTNVQEHPYQGRIIDEK"),p,
  list(enzyme="trypsin",lysine_blocked=FALSE,analysis_scope="specificity_v3"))
stopifnot(all(c("Spanning ends","Spanning terminal form") %in% names(r)),
  match("Spanning ends",names(r))==match("Spanning modifications",names(r))+1L,
  identical(tail(setdiff(names(r),"key"),2),c("Neo form ID","Spanning form ID")))
a <- row_for(r,"TEST1 5 | 6 Neo-N GAPSASVAAAPAPAR")
b <- row_for(r,"TEST1 6 | 7 Neo-N APSASVAAAPAPAR")
c2 <- row_for(r,"TEST2 6 | 7 Neo-N LPDTNVQEHPYQGR")
stopifnot(a$`Spanning sequence`=="FGAPSASVAAAPAPAR",a$`Spanning terminal form`=="No",
  a$`Spanning ends`=="N: digestion-compatible; C: digestion-compatible",
  b$`Spanning sequence`=="FGAPSASVAAAPAPAR",b$`Spanning terminal form`=="No",
  c2$`Spanning sequence`=="GLPDTNVQEHPYQGR",c2$`Spanning terminal form`=="Yes",
  c2$`Spanning ends`=="N: not digestion-compatible; C: digestion-compatible")

# Pre-digestion N-terminal label with lysine blocking: a labelled internal
# N terminus is a terminal form even when it follows R.
seqs <- c("QEDLPETFDAR","EDLPETFDAR","DLPETFDAR","MSTR")
p <- data.frame("Peptide sequence"=seqs,"Modification detail"=c(rep("1xTMT [N-Term]",3),"1xTMT [N-Term]"),
  "Representative ID"="TEST1","Reported IDs"="TEST1",check.names=FALSE)
p <- cbind(p,quant(c(1000,100,100,500),c(250,400,400,500)))
r <- run_case(c(">sp|TEST1|ONE GN=ONE","MSTRQEDLPETFDARLLSER"),p,
  list(enzyme="trypsin",lysine_blocked=TRUE,pre_digest_n_label="TMT",analysis_scope="specificity_v3"))
d <- row_for(r,"TEST1 6 | 7 Neo-N DLPETFDAR")
stopifnot(d$`Spanning terminal form`=="Yes",grepl("^N: pre-digestion N label; C: digestion-compatible$",d$`Spanning ends`))
cat("Spanning ends: terminal spanning forms are reported in Results without changing selection.\n")
