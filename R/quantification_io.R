# Resolve one complete export layout; never combine tables from different runs.
pm_quantification_input_paths <- function(directory, require_protein = FALSE) {
  directory <- normalizePath(directory, winslash="/", mustWork=TRUE)
  required <- c("peptides", "matches", "measurements")
  files <- c(peptides="peptides.csv", matches="matches.csv", measurements="measurements.csv", proteins="source_proteins.csv")
  roots <- c(directory, file.path(directory,"audit"))
  present <- vapply(roots,function(path) any(file.exists(file.path(path,files))),logical(1))
  if(sum(present)>1L) stop("Ambiguous mapping exports: tables occur in both the result directory and audit/. Choose one complete export directory.")
  if(!any(present)) stop("Mapping results require peptides.csv, matches.csv and measurements.csv in the result directory or audit/.")
  selected <- roots[present]
  paths <- setNames(file.path(selected,files),names(files))
  needed <- c(required,if(require_protein)"proteins")
  missing <- needed[!file.exists(paths[needed]) | dir.exists(paths[needed])]
  if(length(missing)) stop("Missing mapping export: ",paste(files[missing],collapse=", "),". Restore the complete export; tables are not combined across directories.")
  paths <- paths[file.exists(paths)]
  actual <- normalizePath(paths,winslash="/",mustWork=TRUE)
  if(any(dir.exists(paths)) || any(!startsWith(actual,paste0(directory,"/"))))
    stop("Mapping export path escapes its result directory.")
  setNames(actual,names(paths))
}

# Read the documented mapping export, without rerunning sequence analysis.
pm_load_quantification_result <- function(directory) {
  paths <- pm_quantification_input_paths(directory)
  required <- c("peptides", "matches", "measurements")
  out <- setNames(lapply(paths[required], pm_read_table), required)
  number <- function(data, columns) {
    for (name in intersect(columns, names(data))) {
      original <- data[[name]]
      value <- suppressWarnings(as.numeric(original))
      if (any(!is.na(original) & nzchar(original) & is.na(value))) stop("Invalid numeric field in mapping export: ", name)
      data[[name]] <- value
    }
    data
  }
  out$peptides <- number(out$peptides, c("source_row", "match_count", "matched_accession_count"))
  out$matches <- number(out$matches, c("record_id", "start", "end", "protein_length"))
  out$measurements <- number(out$measurements, "value")
  for (table in names(out)) for (column in intersect(c("is_master", "included_in_analysis", "modifications_available"), names(out[[table]])))
    out[[table]][[column]] <- pm_flag(out[[table]][[column]])
  if ("proteins" %in% names(paths)) out$source_proteins <- pm_read_table(paths[["proteins"]])
  out
}

pm_export_quantification <- function(quantification, metadata, output_dir, inputs, assumptions = character(), protein_abundance_source = "Not supplied") {
  if (!requireNamespace("openxlsx", quietly=TRUE) || !requireNamespace("jsonlite",quietly=TRUE) || !requireNamespace("digest",quietly=TRUE))
    stop("Install openxlsx, jsonlite and digest for quantitative export.")
  if (dir.exists(output_dir) && length(list.files(output_dir,all.files=TRUE,no..=TRUE))) stop("Choose an empty output directory.")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  tables <- list("Peptide effects"=quantification$effects,"Measurement pairs"=quantification$measurement_pairs,
                 "Samples"=metadata$samples,"Comparisons"=metadata$comparisons)
  for (name in c("effects","measurement_pairs")) pm_write_csv(quantification[[name]],file.path(output_dir,paste0(name,".csv")))
  for (name in names(metadata)) pm_write_csv(metadata[[name]],file.path(output_dir,paste0(name,".csv")))
  priority <- c("comparison","sequence","modifications","accession","raw_peptide_log2fc","matched_peptide_log2fc",
    "protein_log2fc","adjusted_log2fc","status","candidate_accession_count","n_paired_treatment","n_paired_control")
  tables[["Peptide effects"]] <- tables[["Peptide effects"]][,c(intersect(priority,names(quantification$effects)),setdiff(names(quantification$effects),priority)),drop=FALSE]
  labels <- c(comparison="Comparison",sequence="Peptide sequence",modifications="Modifications",accession="PD master accession",
    raw_peptide_log2fc="Raw peptide log2FC",matched_peptide_log2fc="Matched peptide log2FC",protein_log2fc="Protein log2FC",
    adjusted_log2fc="Adjusted log2FC",status="Result",candidate_accession_count="PD master candidates",
    n_paired_treatment="Matched treatment samples",n_paired_control="Matched control samples",
    peptide_id="Peptide ID",source_row="Input row",protein_source_row="Protein input row",
    protein_identified_peptide_count="Protein identified peptides",protein_psm_count="Protein PSMs",protein_quantifying_peptide_count="Protein quantifying peptides")
  workbook_columns <- list()
  for(name in names(tables)) {
    original <- names(tables[[name]])
    display <- original
    if(name %in% c("Peptide effects","Measurement pairs")) {
      display <- gsub("_"," ",original,fixed=TRUE)
      substring(display,1,1) <- chartr("abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ",substring(display,1,1))
      known <- original %in% names(labels);display[known] <- unname(labels[original[known]])
    }
    workbook_columns[[name]] <- as.list(setNames(display,original))
    names(tables[[name]]) <- display
  }
  pm_with_excel_utf8({
    wb <- openxlsx::createWorkbook(creator="Atelente",title="Descriptive protein-adjusted peptide effects")
    header <- openxlsx::createStyle(fontColour="#FFFFFF",fgFill="#263F50",textDecoration="bold",wrapText=TRUE,valign="center")
    openxlsx::modifyBaseFont(wb,fontName="Arial",fontSize=10,fontColour="#172A39")
    for (name in names(tables)) {
      data <- tables[[name]]
      if (nrow(data)>1048575L || ncol(data)>16384L) stop("Excel limits exceeded for ",name,". Complete CSV tables were written.")
      openxlsx::addWorksheet(wb,name,gridLines=FALSE)
      openxlsx::writeData(wb,name,data,headerStyle=header,withFilter=nrow(data)>0,keepNA=FALSE)
      openxlsx::freezePane(wb,name,firstRow=TRUE,firstCol=TRUE)
      openxlsx::setRowHeights(wb,name,1,42)
      for (j in seq_along(data)) {
        is_text <- is.character(data[[j]])
        width <- if(grepl("sequence|modification|column|status|reason",names(data)[j],ignore.case=TRUE)) 38 else 24
        openxlsx::setColWidths(wb,name,j,width)
        format <- if(is_text) "@" else if(grepl("log2FC",names(data)[j],fixed=TRUE)) "0.0000" else "GENERAL"
        if(nrow(data)) openxlsx::addStyle(wb,name,openxlsx::createStyle(numFmt=format,valign="center"),
          rows=seq_len(nrow(data))+1L,cols=j,gridExpand=TRUE)
      }
    }
    openxlsx::addWorksheet(wb,"Analysis notes",gridLines=FALSE)
    notes <- c("Descriptive protein-adjusted peptide effects",
      "Treatment minus Control on the log2 scale. Adjusted effect = matched peptide effect minus matched protein effect.",
      "Peptide and protein effects for adjustment use the same positive, observed samples within each condition. Raw peptide effect is shown separately.",
      "Alternative PD master accessions remain separate candidate estimates. These are not additional independent observations.",
      "No p-values, imputation, site aggregation or spanning-specific correction policy is applied.",
      paste0("Protein abundance source: ",protein_abundance_source),paste0("Assumption: ",as.character(assumptions)))
    openxlsx::writeData(wb,"Analysis notes",data.frame(Notes=notes),colNames=FALSE)
    openxlsx::setColWidths(wb,"Analysis notes",1,110)
    openxlsx::setRowHeights(wb,"Analysis notes",seq_along(notes),42)
    openxlsx::addStyle(wb,"Analysis notes",openxlsx::createStyle(wrapText=TRUE,valign="center"),rows=seq_along(notes),cols=1,gridExpand=TRUE)
    openxlsx::addStyle(wb,"Analysis notes",header,rows=1,cols=1)
    openxlsx::saveWorkbook(wb,file.path(output_dir,"protein_adjusted_effects.xlsx"),overwrite=FALSE)
  })
  manifest <- list(version="0.5.0",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    settings=quantification$settings,protein_abundance_source=protein_abundance_source,assumptions=assumptions,
    workbook_columns=workbook_columns,
    inputs=lapply(inputs,function(path) list(file=normalizePath(path,winslash="/",mustWork=TRUE),sha256=digest::digest(file=path,algo="sha256"))),
    effect_rows=nrow(quantification$effects),measurement_pair_rows=nrow(quantification$measurement_pairs))
  jsonlite::write_json(manifest,file.path(output_dir,"quantification_manifest.json"),pretty=TRUE,auto_unbox=TRUE,na="null",null="null")
  invisible(manifest)
}
