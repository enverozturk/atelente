# Output for the neutral quantity workflow. Quantitative source grain is retained.
pm_export_sample_ratios <- function(data, result, output_dir, inputs, excel=TRUE, notes=NULL) {
  if(!requireNamespace("jsonlite",quietly=TRUE) || !requireNamespace("digest",quietly=TRUE)) stop("Install jsonlite and digest for export.")
  if(!is.logical(excel)||length(excel)!=1L||is.na(excel)) stop("excel must be TRUE or FALSE.")
  if(excel&&!requireNamespace("openxlsx",quietly=TRUE)) stop("Install openxlsx or disable Excel export.")
  if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE))) stop("Choose an empty output directory.")
  tables <- c(list(sample_ratios=result$ratios),data[c("observations","entities","links","peptide_measurements","protein_measurements","samples","source_samples","source_peptides","source_proteins")])
  if(any(!vapply(tables,is.data.frame,logical(1)))) stop("Expected canonical tables are missing.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  for(name in names(tables)) pm_write_csv(tables[[name]],file.path(output_dir,paste0(name,".csv")))
  workbook_columns <- list()
  if(excel) {
    visible <- list("Sample ratios"=tables$sample_ratios,"Peptide measurements"=tables$peptide_measurements,
      "Protein measurements"=tables$protein_measurements,"Samples"=tables$samples)
    priority <- c("sequence","modifications","entity_id","accessions","sample_id","condition","normalized_ratio","log2_ratio","ratio_status","ratio_representation","usable_pair")
    visible[[1]] <- visible[[1]][,c(intersect(priority,names(visible[[1]])),setdiff(names(visible[[1]]),priority)),drop=FALSE]
    labels <- c(sequence="Peptide sequence",modifications="Modifications",observation_id="Observation ID",entity_id="Protein / group ID",
      accessions="Member accessions",sample_id="Sample",condition="Condition",normalized_ratio="Peptide / protein",log2_ratio="Log2 ratio",
      ratio_status="Result",ratio_representation="Linear ratio display",usable_pair="Usable pair",input_value="Input quantity",log2_value="Log2 quantity",
      original_value="Original value text",source_column="Source column",source_row="Source row",entity_source_row="Protein source row")
    pm_with_excel_utf8({
      wb <- openxlsx::createWorkbook(creator="Atelente",title="Per-sample peptide/protein ratios")
      openxlsx::modifyBaseFont(wb,fontName="Arial",fontSize=10,fontColour="#172A39")
      header <- openxlsx::createStyle(fontColour="#FFFFFF",fgFill="#263F50",textDecoration="bold",wrapText=TRUE,valign="center")
      for(name in names(visible)) {
        table <- visible[[name]];keys <- names(table)
        display <- gsub("_"," ",keys,fixed=TRUE)
        substring(display,1,1) <- chartr("abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ",substring(display,1,1))
        known <- keys %in% names(labels);display[known] <- unname(labels[keys[known]])
        workbook_columns[[name]] <- as.list(setNames(display,keys));names(table)<-display
        if(nrow(table)>1048575L || ncol(table)>16384L) stop("Excel limits exceeded; complete CSV tables were written. Rerun with excel=false in a new folder.")
        if(any(vapply(table,function(x)is.character(x)&&any(nchar(x,allowNA=TRUE)>32767L,na.rm=TRUE),logical(1)))) stop("Excel cell text limit exceeded; CSV tables were written.")
        openxlsx::addWorksheet(wb,name,gridLines=FALSE)
        openxlsx::writeData(wb,name,table,headerStyle=header,withFilter=nrow(table)>0,keepNA=FALSE)
        openxlsx::freezePane(wb,name,firstRow=TRUE,firstCol=TRUE)
        openxlsx::setRowHeights(wb,name,1,36)
        height <- rep(17,nrow(table))
        for(j in seq_along(table)) {
          long <- grepl("sequence|modification|accession|original|column|reason",keys[j],ignore.case=TRUE)
          width <- if(long) 40 else 24
          openxlsx::setColWidths(wb,name,j,width)
          if(nrow(table)) {
            fmt <- if(is.character(table[[j]])) "@" else "GENERAL"
            openxlsx::addStyle(wb,name,openxlsx::createStyle(numFmt=fmt,wrapText=long,valign="center"),rows=seq_len(nrow(table))+1L,cols=j,gridExpand=TRUE)
            if(long&&is.character(table[[j]])) height <- pmax(height,pmin(409,14*pmax(1,ceiling(nchar(ifelse(is.na(table[[j]]),"",table[[j]]))/(width*.7)))+3))
          }
        }
        if(nrow(table)) openxlsx::setRowHeights(wb,name,seq_len(nrow(table))+1L,height)
      }
      explanation <- c("Per-sample peptide/protein ratios",
        "Each ratio uses the explicitly linked protein or protein-group quantity from the same sample.",
        "Protein groups retain one quantity per sample. Alternative links are separate possible baselines, not independent peptide observations.",
        "Input scales are declared. Log2 zero and negative values are valid. Missing, zero and negative linear inputs are not imputed.",
        "Finite log2 ratios can remain usable when a linear ratio is too large or too small to represent; the display field explains this.",
        "No condition comparison, p-value, cleavage fraction or protease activity is inferred. FASTA placement is not performed by this command.",
        paste0("Peptide quantity: ",data$spec$peptide$quantity_type,"; scale: ",data$spec$peptide$scale),
        paste0("Protein quantity: ",data$spec$protein$quantity_type,"; scale: ",data$spec$protein$scale),as.character(unlist(notes)))
      openxlsx::addWorksheet(wb,"Read me",gridLines=FALSE)
      openxlsx::writeData(wb,"Read me",data.frame(Notes=explanation),colNames=FALSE)
      openxlsx::setColWidths(wb,"Read me",1,110)
      openxlsx::setRowHeights(wb,"Read me",seq_along(explanation),42)
      openxlsx::addStyle(wb,"Read me",openxlsx::createStyle(wrapText=TRUE,valign="center"),rows=seq_along(explanation),cols=1,gridExpand=TRUE)
      openxlsx::addStyle(wb,"Read me",header,rows=1,cols=1)
      openxlsx::saveWorkbook(wb,file.path(output_dir,"sample_ratios.xlsx"),overwrite=FALSE)
    })
  }
  manifest <- list(version="0.6.1",stage="per_sample_ratios",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    spec=data$spec,settings=result$settings,notes=notes,excel_export=excel,workbook_columns=workbook_columns,
    row_counts=as.list(vapply(tables,nrow,integer(1))),
    inputs=lapply(inputs,function(path)list(file=normalizePath(path,winslash="/",mustWork=TRUE),sha256=digest::digest(file=path,algo="sha256"))),
    csv_sha256=as.list(setNames(vapply(names(tables),function(name)digest::digest(file=file.path(output_dir,paste0(name,".csv")),algo="sha256"),character(1)),names(tables))))
  jsonlite::write_json(manifest,file.path(output_dir,"normalization_manifest.json"),auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
  invisible(manifest)
}
