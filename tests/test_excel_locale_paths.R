# All workbook-producing entry points must share the character-locale guard.
# Fixture literals use escapes so LC_ALL=C can parse the test itself.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
for(module in c("report","excel_locale","sample_metadata","quantification_io","normalization_io","statistics_io"))
  source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")

check_paths <- function() {
  original <- Sys.getlocale("LC_CTYPE")
  on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE",original)),add=TRUE)
  numeric <- Sys.getlocale("LC_NUMERIC")
  td <- tempfile("excel-locale-paths-");dir.create(td)
  label <- "\u0130la\u00e7 \u03b1 & <1>"
  heading <- "\u96f6 \U0001F9EA"
  values <- data.frame(text=c(label,NA_character_,"00123"),quantity=c(pi,0,NA_real_),check.names=FALSE)
  names(values)[1] <- heading
  samples <- data.frame(Sample=label,Condition=label,Use="Yes",
    "Peptide column"=label,"Protein column"=heading,check.names=FALSE)
  metadata <- list(samples=samples,comparisons=data.frame(Name=label,Treatment=label,Control="Control"))
  quant <- list(effects=values,measurement_pairs=values,settings=list())
  data <- setNames(rep(list(values),9),c("observations","entities","links","peptide_measurements",
    "protein_measurements","samples","source_samples","source_peptides","source_proteins"))
  data$spec <- list(peptide=list(quantity_type=label,scale="linear"),protein=list(quantity_type=label,scale="linear"))
  ratios <- list(ratios=values,settings=list())
  statistics <- list(results=values,measurements=values,samples=samples,comparisons=metadata$comparisons,settings=list())
  snapshot <- serialize(list(metadata,quant,data,ratios,statistics),NULL)
  input <- file.path(td,"input.txt");writeLines("synthetic input",input)
  writers <- list(
    metadata=function(dir) {dir.create(dir);pm_write_sample_metadata(metadata,file.path(dir,"metadata.xlsx"))},
    quantification=function(dir) pm_export_quantification(quant,metadata,dir,list(input=input),assumptions=label),
    normalization=function(dir) pm_export_sample_ratios(data,ratios,dir,list(input=input),notes=label),
    statistics=function(dir) pm_export_statistics(statistics,dir,list(input=input),notes=label))
  cells <- function(path) lapply(setNames(readxl::excel_sheets(path),readxl::excel_sheets(path)),
    function(sheet) as.data.frame(readxl::read_xlsx(path,sheet=sheet,.name_repair="minimal")))
  for(name in names(writers)) {
    stopifnot(nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE","C"))))
    c_dir <- file.path(td,paste0(name,"-c"));writers[[name]](c_dir)
    stopifnot(Sys.getlocale("LC_CTYPE")=="C",Sys.getlocale("LC_NUMERIC")==numeric)
    pm_excel_utf8_ctype();utf8 <- Sys.getlocale("LC_CTYPE")
    u_dir <- file.path(td,paste0(name,"-utf8"));writers[[name]](u_dir)
    stopifnot(Sys.getlocale("LC_CTYPE")==utf8,Sys.getlocale("LC_NUMERIC")==numeric)
    a <- cells(list.files(c_dir,pattern="\\.xlsx$",full.names=TRUE))
    b <- cells(list.files(u_dir,pattern="\\.xlsx$",full.names=TRUE))
    stopifnot(identical(a,b),any(vapply(a,function(x) any(unlist(x)==label,na.rm=TRUE),logical(1))))
    if(name!="metadata") stopifnot(heading %in% names(a[[1]]),
      identical(a[[1]][[heading]],values[[heading]]),
      isTRUE(all.equal(a[[1]][[2]],values$quantity,tolerance=1e-14)))
    csv <- list.files(c_dir,pattern="\\.csv$")
    stopifnot(identical(csv,list.files(u_dir,pattern="\\.csv$")))
    for(file in csv) stopifnot(identical(digest::digest(file=file.path(c_dir,file),algo="sha256"),
                                       digest::digest(file=file.path(u_dir,file),algo="sha256")))
    # Refusal must not create a workbook, and must restore the caller locale.
    candidates <- pm_excel_utf8_candidates
    assign("pm_excel_utf8_candidates",character(),envir=environment(pm_excel_utf8_ctype))
    stopifnot(nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE","C"))))
    refused_dir <- file.path(td,paste0(name,"-refused"))
    err <- try(writers[[name]](refused_dir),silent=TRUE)
    assign("pm_excel_utf8_candidates",candidates,envir=environment(pm_excel_utf8_ctype))
    stopifnot(inherits(err,"try-error"),grepl("UTF-8 locale",err,fixed=TRUE),
      !length(list.files(refused_dir,pattern="\\.xlsx$")),Sys.getlocale("LC_CTYPE")=="C")
  }
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",original)))
  stopifnot(identical(snapshot,serialize(list(metadata,quant,data,ratios,statistics),NULL)))
  cat("All four workbook paths: Unicode cells/headers and numeric cells preserved; CSV bytes identical; locale restored; no workbook on refusal.\n")
}
check_paths()
