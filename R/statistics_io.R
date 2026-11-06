# Load a verified ratio export; identifiers remain literal strings.
pm_load_ratio_statistics <- function(directory) {
  manifest_path <- file.path(directory,"normalization_manifest.json")
  ratio_path <- file.path(directory,"sample_ratios.csv")
  if(!file.exists(manifest_path)||!file.exists(ratio_path)) stop("Require a normalization manifest and sample_ratios.csv.")
  manifest <- jsonlite::read_json(manifest_path,simplifyVector=TRUE)
  if(!identical(manifest$stage,"per_sample_ratios")) stop("Expected per_sample_ratios input stage.")
  expected <- manifest$csv_sha256$sample_ratios
  if(is.null(expected)||!identical(expected,digest::digest(file=ratio_path,algo="sha256"))) stop("Ratio CSV hash differs from its manifest.")
  ratios <- pm_read_table(ratio_path)
  # CSV blanks are imported as NA; the canonical unknown condition is an empty string.
  if("condition" %in% names(ratios)) ratios$condition[is.na(ratios$condition)] <- ""
  number <- c("source_row","entity_source_row","candidate_entity_count","peptide_input_value","protein_input_value",
    "peptide_log2_value","protein_log2_value","normalized_ratio","log2_ratio")
  boolean <- c("modifications_available","observation_excluded","entity_excluded","usable_pair")
  for(name in intersect(number,names(ratios))) {
    original <- ratios[[name]];value <- suppressWarnings(as.numeric(original))
    if(any(!is.na(original)&!is.finite(value))) stop("Invalid numeric ratio field: ",name)
    ratios[[name]] <- value
  }
  for(name in intersect(boolean,names(ratios))) {
    original <- ratios[[name]]
    if(any(!is.na(original)&!original %in% c("TRUE","FALSE"))) stop("Invalid logical ratio field: ",name)
    ratios[[name]] <- ifelse(is.na(original),NA,original=="TRUE")
  }
  list(ratios=ratios,manifest=manifest,inputs=list(ratio_manifest=manifest_path,ratios=ratio_path))
}

pm_export_statistics <- function(result,output_dir,inputs,excel=TRUE,notes=NULL) {
  if(!is.logical(excel)||length(excel)!=1L||is.na(excel)) stop("excel must be TRUE or FALSE.")
  if(excel&&!requireNamespace("openxlsx",quietly=TRUE)) stop("Install openxlsx or disable Excel export.")
  if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE))) stop("Choose an empty output directory.")
  tables <- result[c("results","measurements","samples","comparisons")]
  if(length(tables)!=4L||any(!vapply(tables,is.data.frame,logical(1)))) stop("Expected statistics tables are missing.")
  for(table in tables) for(column in table) if(is.numeric(column)&&any(!is.na(column)&!is.finite(column))) stop("Encode infinite values explicitly before export.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  for(name in names(tables)) pm_write_csv(tables[[name]],file.path(output_dir,paste0(name,".csv")))
  workbook_columns <- list()
  if(excel) {
    visible <- setNames(tables,c("Results","Sample values","Samples","Comparisons"))
    priority <- c("comparison","sequence","modifications","entity_id","accessions","log2_effect","p_value","adjusted_p_value","status","n_control","n_treatment")
    visible[[1]] <- visible[[1]][,c(intersect(priority,names(visible[[1]])),setdiff(names(visible[[1]]),priority)),drop=FALSE]
    labels <- c(comparison="Comparison",sequence="Peptide sequence",modifications="Modifications",entity_id="Protein / group ID",
      observation_id="Observation ID",accessions="Member accessions",log2_effect="Adjusted log2 effect",p_value="P-value",
      adjusted_p_value="BH-adjusted P-value",status="Result",n_control="Usable control units",n_treatment="Usable treatment units",
      log2_ratio="Log2 peptide / protein",prior_df_infinite="Infinite prior degrees of freedom")
    pm_with_excel_utf8({
      wb <- openxlsx::createWorkbook(creator="Atelente",title="Protein-adjusted peptide statistics")
      openxlsx::modifyBaseFont(wb,fontName="Arial",fontSize=10,fontColour="#172A39")
      header <- openxlsx::createStyle(fontColour="#FFFFFF",fgFill="#263F50",textDecoration="bold",wrapText=TRUE,valign="center")
      for(name in names(visible)) {
        table <- visible[[name]];keys <- names(table)
        display <- gsub("_"," ",keys,fixed=TRUE)
        substring(display,1,1) <- chartr("abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ",substring(display,1,1))
        known <- keys %in% names(labels);display[known] <- unname(labels[keys[known]])
        workbook_columns[[name]] <- as.list(setNames(display,keys));names(table) <- display
        if(nrow(table)>1048575L||ncol(table)>16384L) stop("Excel limits exceeded; full CSV tables were written.")
        if(any(vapply(table,function(x)is.character(x)&&any(nchar(x)>32767,na.rm=TRUE),logical(1)))) stop("Excel cell limit exceeded; CSV tables were written.")
        openxlsx::addWorksheet(wb,name,gridLines=FALSE)
        openxlsx::writeData(wb,name,table,headerStyle=header,withFilter=nrow(table)>0,keepNA=FALSE)
        openxlsx::freezePane(wb,name,firstRow=TRUE,firstCol=TRUE);openxlsx::setRowHeights(wb,name,1,44)
        height <- rep(18,nrow(table))
        for(j in seq_along(table)) {
          long <- grepl("sequence|modification|accession|reason|status|column",keys[j],ignore.case=TRUE)
          width <- if(long) 40 else 26
          openxlsx::setColWidths(wb,name,j,width)
          if(nrow(table)) {
            format <- if(is.character(table[[j]])) "@" else "GENERAL"
            openxlsx::addStyle(wb,name,openxlsx::createStyle(numFmt=format,wrapText=long,valign="center"),rows=seq_len(nrow(table))+1L,cols=j,gridExpand=TRUE)
            if(long&&is.character(table[[j]])) height <- pmax(height,pmin(409,14*pmax(1,ceiling(nchar(ifelse(is.na(table[[j]]),"",table[[j]]))/(width*.7)))+3))
          }
        }
        if(nrow(table)) openxlsx::setRowHeights(wb,name,seq_len(nrow(table))+1L,height)
      }
      explanation <- c("Protein-adjusted peptide statistics",
        "Response: log2(peptide quantity / protein quantity). Positive effects mean higher in Treatment than Control.",
        "Independent experimental units only; repeated measurements of one unit are not independent replicates. Paired and batch designs are not fitted.",
        "Each contrast fits a two-condition limma model with pooled within-condition variance, trend=false and robust=false.",
        "BH adjustment is within each contrast across tested observations. It does not control the combined family across several contrasts.",
        "Only observations with one declared protein/group candidate enter inference. Alternative denominators remain diagnostic rows.",
        "At least two usable units per condition and two testable observations are required. Exactly zero or invalid residual variances are not tested.",
        "There is no imputation or significance cutoff. Missingness related to abundance or condition may bias observed-value inference.",
        "Infinite prior degrees of freedom is a valid limma result: its numeric cell is blank and the explicit infinite flag is true.",
        "Shared protein baselines and quantities derived from the same peptides create dependencies. These results are not site-level cleavage or protease-activity tests.",
        as.character(unlist(notes)))
      openxlsx::addWorksheet(wb,"Read me",gridLines=FALSE)
      openxlsx::writeData(wb,"Read me",data.frame(Notes=explanation),colNames=FALSE)
      openxlsx::setColWidths(wb,"Read me",1,110);openxlsx::setRowHeights(wb,"Read me",seq_along(explanation),44)
      openxlsx::addStyle(wb,"Read me",openxlsx::createStyle(wrapText=TRUE,valign="center"),rows=seq_along(explanation),cols=1,gridExpand=TRUE)
      openxlsx::addStyle(wb,"Read me",header,rows=1,cols=1)
      openxlsx::saveWorkbook(wb,file.path(output_dir,"peptide_statistics.xlsx"),overwrite=FALSE)
    })
  }
  manifest <- list(version="0.7.0",stage="independent_unit_statistics",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    settings=result$settings,notes=notes,excel_export=excel,workbook_columns=workbook_columns,
    row_counts=as.list(vapply(tables,nrow,integer(1))),
    inputs=lapply(inputs,function(path)list(file=normalizePath(path,winslash="/",mustWork=TRUE),sha256=digest::digest(file=path,algo="sha256"))),
    csv_sha256=as.list(setNames(vapply(names(tables),function(name)digest::digest(file=file.path(output_dir,paste0(name,".csv")),algo="sha256"),character(1)),names(tables))))
  jsonlite::write_json(manifest,file.path(output_dir,"statistics_manifest.json"),auto_unbox=TRUE,pretty=TRUE,na="null",null="null",digits=I(17))
  invisible(manifest)
}
