# Reader-facing typed Excel views. Original scientific tables remain unchanged.
pm_excel_reference <- function(x) {
  vapply(as.character(x), function(value) {
    if (is.na(value) || !nzchar(value)) return(value)
    # Only the inspected complete MEROPS wrapper is formatted. Unknown syntax
    # remains verbatim; reference IDs and raw source fields remain in CSV.
    match <- regmatches(value, regexec("^<%([^<>]*)\\[([A-Za-z0-9]+)\\]%>$", value, perl=TRUE))[[1]]
    if (!length(match)) return(value)
    trimws(gsub("%etal", "et al.", match[2], fixed=TRUE))
  }, character(1), USE.NAMES=FALSE)
}

pm_excel_select <- function(data, fields) {
  if (is.null(data)) return(NULL)
  keep <- intersect(names(fields), names(data))
  out <- data[,keep,drop=FALSE]
  names(out) <- unname(fields[keep])
  out
}

pm_site_summary <- function(result) {
  sites <- result$cleavage$sites
  if (is.null(sites)) return(data.frame(site_id=character(),record_id=integer(),accession=character(),
    bond_left=integer(),bond_right=integer(),context_display=character(),context_complete=logical(),
    neo_n_count=integer(),neo_c_count=integer(),spanning_count=integer(),uniprot_summary=character(),
    uniprot_status=character(),merops_summary=character(),merops_reference_count=integer(),merops_motif_count=integer(),
    merops_mismatch_count=integer(),merops_context_unknown_count=integer()))
  fields <- c("site_id","record_id","accession","bond_left","bond_right","context_display","context_complete",
    "neo_n_count","neo_c_count","spanning_count")
  out <- sites[,fields,drop=FALSE]
  out$uniprot_summary <- rep("Not requested",nrow(out));out$uniprot_status <- rep("not_requested",nrow(out))
  out$merops_summary <- rep("Not requested",nrow(out));out$merops_reference_count<-rep(NA_integer_,nrow(out));out$merops_motif_count<-rep(NA_integer_,nrow(out))
  out$merops_mismatch_count<-rep(NA_integer_,nrow(out));out$merops_context_unknown_count<-rep(NA_integer_,nrow(out))
  if (!is.null(result$uniprot)) {
    ur<-result$uniprot$records;ua<-result$uniprot$site_annotations
    index<-split(seq_len(nrow(ua)),ua$site_id)
    for(i in seq_len(nrow(out))) {
      status<-ur$status[match(out$record_id[i],ur$record_id)]
      if(length(status)!=1L || is.na(status)) status<-"not_assessed"
      out$uniprot_status[i]<-status
      labels<-c(matched="No matching boundary annotation",not_found="No entry in this snapshot",
        unsupported_identifier="Not queried: non-UniProt identifier",
        sequence_mismatch="Not applied: reference sequence differs",sequence_unavailable="Not applied: reference sequence unavailable",
        unsupported_isoform="Not applied: isoform reference unavailable",ambiguous_accession="Not applied: ambiguous reference",
        not_assessed="Reference not assessed")
      out$uniprot_summary[i]<-if(status %in% names(labels)) unname(labels[status]) else status
      rows<-index[[out$site_id[i]]]
      if(length(rows)) out$uniprot_summary[i]<-paste(unique(paste(ua$feature_type[rows],gsub("_"," ",ua$boundary[rows],fixed=TRUE),sep=": ")),collapse="; ")
    }
  }
  if(!is.null(result$merops)) {
    ma<-result$merops$matches;index<-split(seq_len(nrow(ma)),ma$site_id)
    for(i in seq_len(nrow(out))) {
      rows<-index[[out$site_id[i]]]
      a<-sum(ma$relation_type[rows]=="reference_site_candidate");b<-sum(ma$relation_type[rows]=="motif_compatibility")
      coordinate_rows<-rows[ma$relation_type[rows]=="reference_site_candidate"]
      mismatch<-sum(ma$context_status[coordinate_rows]=="mismatch");unknown<-sum(ma$context_status[coordinate_rows]=="unknown")
      out$merops_reference_count[i]<-as.integer(a);out$merops_motif_count[i]<-as.integer(b)
      out$merops_mismatch_count[i]<-as.integer(mismatch);out$merops_context_unknown_count[i]<-as.integer(unknown)
      qualifiers<-c(if(mismatch) paste(mismatch,"context mismatch"),if(unknown) paste(unknown,"context unknown"))
      out$merops_summary[i]<-if(!length(rows)) "No reference records" else paste0(a," accession/position records",
        if(length(qualifiers)) paste0(" (",paste(qualifiers,collapse=", "),")"),"; ",b," motif-only records")
    }
  }
  out
}

pm_excel_tables <- function(result) {
  tables<-list()
  if(!is.null(result$cleavage)) tables[["Cleavage sites"]]<-pm_excel_select(pm_site_summary(result),c(
    accession="Protein accession",bond_left="Bond left",bond_right="Bond right",context_display="P4-P4 prime",
    neo_n_count="Neo-N candidates",neo_c_count="Neo-C candidates",spanning_count="Spanning observations",
    uniprot_summary="UniProt processing",merops_summary="MEROPS reference summary",context_complete="Complete context",
    uniprot_status="UniProt reference status",merops_reference_count="Accession-position records",merops_motif_count="Motif-only records",
    merops_mismatch_count="Context mismatch records",merops_context_unknown_count="Context unknown records",
    record_id="FASTA record",site_id="Site ID"))
  if(!is.null(result$cleavage$evidence)) {
    data<-result$cleavage$evidence
    if(anyDuplicated(result$peptides$peptide_id)) stop("Peptide IDs must be unique for Excel evidence lookup.")
    obs<-match(data$peptide_id,result$peptides$peptide_id)
    for(field in c("sequence","modifications","source_row")) data[[field]]<-result$peptides[[field]][obs]
    tables[["Cleavage evidence"]]<-pm_excel_select(data,c(accession="Protein accession",bond_left="Bond left",
      bond_right="Bond right",classification="Classification",role="Peptide terminus",sequence="Peptide sequence",
      modifications="Modifications",start="Peptide start",end="Peptide end",source_row="Input row",peptide_id="Peptide ID",
      record_id="FASTA record",placement_id="Placement ID",site_id="Site ID",candidate_reason="Candidate reason",
      digestion_compatible="Digestion compatible",pre_digest_n_label="Pre-digest N label",is_candidate="Candidate"))
  }
  tables[["Peptides"]]<-pm_excel_select(result$peptides,c(sequence="Peptide sequence",modifications="Modifications",
    mapping_status="Mapping result",match_count="Position count",matched_accession_count="Protein accession count",
    source_row="Input row",peptide_id="Peptide ID",annotated_sequence="Original sequence annotation",
    modifications_available="Modifications column supplied",exclusion_reason="Exclusion reason"))
  tables[["Protein positions"]]<-pm_excel_select(result$matches,c(accession="Protein accession",start="Start",end="End",
    peptide_id="Peptide ID",record_id="FASTA record",protein_length="Protein length",is_master="PD master protein",
    pd_coordinate_agrees="Matches PD position",reported_by_pd="Reported by PD",left_flank="Left flank",right_flank="Right flank"))
  tables[["Abundances"]]<-pm_excel_select(result$measurements,c(peptide_id="Peptide ID",measurement_label="Measurement label",
    value="Abundance",value_status="Value status",abundance_type="Abundance type",source_column="Source column",
    original_value="Original value text",included_in_analysis="Included in analysis"))
  if(!is.null(result$merops)) {
    data<-result$merops$matches
    data$readable_reference<-pm_excel_reference(data$reference_raw)
    data$readable_relation<-ifelse(data$relation_type=="reference_site_candidate","Accession + position",
      ifelse(data$relation_type=="motif_compatibility","Motif only",data$relation_type))
    tables[["MEROPS annotations"]]<-pm_excel_select(data,c(site_accession="Protein accession",site_bond_left="Bond left",
      site_bond_right="Bond right",site_context="P4-P4 prime",protease_id="Protease ID",protease_name="Protease name",
      readable_relation="Reference relation",context_status="Context agreement",cleavage_type="Source evidence type",
      readable_reference="Reference",reference_ids="Reference IDs",merops_record_id="MEROPS record ID",
      substrate_accession="Reference substrate accession",p1_position="Reference P1 position text",substrate_species="Substrate species",
      method_code_raw="Source method code",additional_evidence_code_raw="Source additional evidence code",site_id="Site ID"))
  }
  if(!is.null(result$uniprot)) {
    fields<-c(accession="Protein accession",bond_left="Bond left",bond_right="Bond right",feature_type="Feature",
      boundary="Boundary",description="Description",evidence_codes="Evidence codes",uniprot_accession="UniProt accession",
      source_feature_id="UniProt feature ID",record_id="FASTA record",feature_key="Feature key",site_id="Site ID")
    tables[["UniProt sites"]]<-pm_excel_select(result$uniprot$site_annotations,fields)
    tables[["UniProt termini"]]<-pm_excel_select(result$uniprot$terminal_annotations,c(
      peptide_id="Peptide ID",accession="Protein accession",role="Peptide terminus",start="Peptide start",end="Peptide end",
      fields[!names(fields) %in% c("accession","site_id")],placement_id="Placement ID"))
    tables[["Reference checks"]]<-pm_excel_select(result$uniprot$records,c(accession="Protein accession",status="UniProt reference status",
      uniprot_accession="UniProt accession",entry_version="Entry version text",sequence_version="Sequence version text",
      annotation_date="Annotation date text",sequence_date="Sequence date text",record_id="FASTA record",
      fasta_sequence_sha256="FASTA sequence SHA256",uniprot_sequence_sha256="UniProt sequence SHA256"))
  }
  # Analysis tables retain their native row grain; alternative mapping rows are
  # distinct from the deduplicated measurement-family test table.
  if (!is.null(result$analysis_tables)) {
    labels <- c(neo_spanning_contrasts="Pair contrasts",pair_results="Pair mapping results",
      pair_measurement_results="Pair tests",pair_excluded="Pair exclusions",
      placement_specificity="Specificity",excluded_placements="Excluded placements")
    for (name in intersect(names(labels),names(result$analysis_tables)))
      tables[[labels[[name]]]] <- result$analysis_tables[[name]]
  }
  # Only valid presentation tables. Each remains at its source grain.
  tables<-tables[!vapply(tables,is.null,logical(1))]
  for(name in names(tables)) {
    data<-tables[[name]]
    for(column in names(data)) {
      if(is.factor(data[[column]])) data[[column]]<-as.character(data[[column]])
      if(grepl(" ID$| IDs$| accession$|^FASTA record$| text$",column)) data[[column]]<-as.character(data[[column]])
      if(is.list(data[[column]]) || !is.atomic(data[[column]])) stop("Unsupported Excel column: ",name," / ",column)
    }
    if(nrow(data)>1048575L || ncol(data)>16384L) stop("Excel worksheet limits exceeded for ",name,"; use the complete CSV tables.")
    if(any(vapply(data,function(x) is.character(x) && any(nchar(x,type="chars",allowNA=TRUE)>32767L,na.rm=TRUE),logical(1))))
      stop("Excel cell text limit exceeded in ",name,"; use the complete CSV tables.")
    tables[[name]]<-data
  }
  tables
}

pm_write_excel <- function(result,path) {
  pm_write_excel_tables(pm_excel_tables(result),path)
}

pm_write_excel_tables <- function(tables,path,sheets=names(tables)) {
  if(!requireNamespace("openxlsx",quietly=TRUE)) stop("Install openxlsx for Excel export, or explicitly disable Excel export.")
  if(file.exists(path)) stop("Excel output already exists; choose a new path.")
  if(anyDuplicated(sheets)||any(!sheets %in% names(tables))) stop("Unknown or duplicate workbook sheets.")
  tables<-tables[sheets]
  if(!length(tables)) stop("No tables available for Excel export.")
  # Covers createWorkbook/writeData/saveWorkbook/repair below (R/excel_locale.R).
  previous_ctype<-pm_excel_utf8_ctype()
  if(!is.null(previous_ctype)) on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE",previous_ctype)),add=TRUE)
  wb<-openxlsx::createWorkbook(creator="Atelente",title="Atelente results")
  openxlsx::modifyBaseFont(wb,fontSize=10,fontName="Arial",fontColour="#172A39")
  header<-openxlsx::createStyle(fontColour="#FFFFFF",fgFill="#263F50",textDecoration="bold",halign="center",valign="center",wrapText=TRUE)
  textstyle<-openxlsx::createStyle(numFmt="@",valign="center")
  wrappedtextstyle<-openxlsx::createStyle(numFmt="@",valign="center",wrapText=TRUE)
  numberstyle<-openxlsx::createStyle(numFmt="GENERAL",halign="right",valign="center")
  integerstyle<-openxlsx::createStyle(numFmt="0",halign="right",valign="center")
  # Display precision only: exported numeric values retain their full precision.
  effectstyle<-openxlsx::createStyle(numFmt='[>=0.0005]0.000;[<=-0.0005]-0.000;"0.000"' ,halign="right",valign="center")
  pstyle<-openxlsx::createStyle(numFmt="0.00E+00",halign="right",valign="center")
  for(name in names(tables)) {
    data<-tables[[name]]
    openxlsx::addWorksheet(wb,name,gridLines=FALSE,tabColour=if(name %in% c("Results","Cleavage sites")) "#263F50" else NULL)
    # writeData writes ordinary character values as literal shared strings, not formulas.
    openxlsx::writeData(wb,name,data,colNames=TRUE,rowNames=FALSE,headerStyle=header,withFilter=nrow(data)>0,keepNA=FALSE)
    openxlsx::freezePane(wb,name,firstRow=TRUE,firstCol=TRUE)
    openxlsx::setRowHeights(wb,name,rows=1,heights=34)
    rowlines<-rep(1L,nrow(data));has_wrapped_column<-FALSE
    for(j in seq_along(data)) {
      col<-data[[j]];label<-names(data)[j]
      wrap<-is.character(col) && grepl("summary|processing|^Reference$|Description|sequence$|Modifications$|^(Protein|Gene|Bond|Cleavage window|MEROPS|UniProt|Value)$|^Candidate reason$|^Evidence codes$|^Original sequence annotation$",label,ignore.case=TRUE)
      width<-if(label=="Value")70 else if(label %in% c("MEROPS","UniProt"))45 else
        if(label=="Bond")36 else if(label %in% c("Protein","Gene","Cleavage window"))30 else
        if(grepl("summary|processing|Reference$|Description|sequence$|Modifications|Candidate reason|Evidence codes|SHA256|Site ID|Source column|Original",label,ignore.case=TRUE)) 42 else
        if(is.numeric(col)||is.logical(col)) max(16,min(24,nchar(label)+2)) else max(18,min(30,nchar(label)+2))
      openxlsx::setColWidths(wb,name,cols=j,widths=width)
      if(nrow(data)) {
        style<-if(is.character(col)) { if(wrap) wrappedtextstyle else textstyle } else if(is.numeric(col)) {
          if(grepl("log2FC$",label)) effectstyle else
          if(label %in% c("P-value","BH-adjusted P")) pstyle else
          if(is.integer(col)) integerstyle else numberstyle
        } else NULL
        if(!is.null(style)) openxlsx::addStyle(wb,name,style,rows=seq_len(nrow(data))+1L,cols=j,gridExpand=TRUE,stack=TRUE)
        if(wrap) {
          # Allow spare horizontal room for proportional letters and word wrapping.
          # Heights affect presentation only; technical IDs remain single-line text.
          lines<-vapply(col,function(value) {
            if(is.na(value)||!nzchar(value)) return(1L)
            parts<-strsplit(gsub("\r\n?","\n",value,perl=TRUE),"\n",fixed=TRUE)[[1]]
            as.integer(max(1,sum(pmax(1,ceiling(nchar(parts,type="chars")/(width*0.70))))))
          },integer(1))
          rowlines<-pmax(rowlines,lines);has_wrapped_column<-TRUE
        }
      }
    }
    if(nrow(data)&&has_wrapped_column)
      openxlsx::setRowHeights(wb,name,rows=seq_len(nrow(data))+1L,heights=pmin(409,pmax(15,rowlines*14+3)))
  }
  # Internal catalog backs offline filtered downloads, not the default workbook.
  if("__Results" %in% names(tables)){
    visibility<-openxlsx::sheetVisibility(wb);visibility[names(tables)=="__Results"]<-"hidden"
    openxlsx::sheetVisibility(wb)<-visibility
  }
  openxlsx::saveWorkbook(wb,path,overwrite=FALSE)
  pm_excel_repair_data_parts(path,tables)
  invisible(list(file=basename(path),sheets=names(tables),rows=vapply(tables,nrow,integer(1))))
}

# openxlsx 4.2.5 may emit dimension="A1" and placeholder drawing/VML
# relationships in data-only sheets even though no such parts exist. Repair
# only those generated metadata parts; all cells, types and styles stay intact.
pm_excel_repair_data_parts <- function(path,tables) {
  if(!requireNamespace("zip",quietly=TRUE))stop("Install zip for Excel export.")
  folder<-tempfile("atelente_xlsx_parts_");dir.create(folder)
  on.exit({
    if(dir.exists(folder)) {
      resolved<-normalizePath(folder,winslash="/",mustWork=TRUE)
      temporary_root<-paste0(normalizePath(tempdir(),winslash="/",mustWork=TRUE),"/")
      if(startsWith(resolved,temporary_root)&&startsWith(basename(resolved),"atelente_xlsx_parts_"))
        unlink(resolved,recursive=TRUE,force=TRUE)
    }
  },add=TRUE)
  utils::unzip(path,exdir=folder)
  text<-function(p)paste(readLines(p,warn=FALSE,encoding="UTF-8"),collapse="")
  save<-function(value,p)writeLines(enc2utf8(value),p,useBytes=TRUE)
  attr<-function(tag,name){hit<-regmatches(tag,regexec(paste0('\\b',name,'="([^"]*)"'),tag,perl=TRUE))[[1]];if(length(hit))hit[2] else ""}
  for(i in seq_along(tables)){
    sheet<-file.path(folder,"xl","worksheets",paste0("sheet",i,".xml"));xml<-text(sheet)
    last<-paste0(openxlsx::int2col(max(1L,ncol(tables[[i]]))),nrow(tables[[i]])+1L)
    xml<-sub('<dimension\\b[^>]*/>',paste0('<dimension ref="A1:',last,'"/>'),xml,perl=TRUE)
    save(xml,sheet)
    relationships<-file.path(folder,"xl","worksheets","_rels",paste0("sheet",i,".xml.rels"))
    if(!file.exists(relationships))next
    rel<-text(relationships);tags<-regmatches(rel,gregexpr('<Relationship\\b[^>]*/>',rel,perl=TRUE))[[1]]
    for(tag in tags){
      if(attr(tag,"TargetMode")=="External")next
      target<-file.path(dirname(sheet),attr(tag,"Target"))
      if(file.exists(target))next
      if(!grepl('/(drawing|vmlDrawing)$',attr(tag,"Type")))stop("Missing workbook relationship target: ",attr(tag,"Target"))
      id<-attr(tag,"Id")
      if(grepl(paste0('r:id="',id,'"'),xml,fixed=TRUE))stop("Data workbook unexpectedly references a missing drawing.")
      rel<-sub(tag,"",rel,fixed=TRUE)
    }
    save(rel,relationships)
  }
  files<-list.files(folder,recursive=TRUE,all.files=TRUE,no..=TRUE)
  archive<-tempfile(fileext=".xlsx");on.exit(unlink(archive),add=TRUE)
  zip::zipr(archive,files=files,root=folder,include_directories=FALSE,mode="mirror")
  if(!file.copy(archive,path,overwrite=TRUE))stop("Could not finalize Excel metadata.")
  invisible(path)
}
