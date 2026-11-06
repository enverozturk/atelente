args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1]),winslash="/")))
source(file.path(root,"R","excel_locale.R"),encoding="UTF-8");source(file.path(root,"R","excel.R"),encoding="UTF-8")
stopifnot(requireNamespace("openxlsx",quietly=TRUE))

# Small independent fixture: multiple annotations must not multiply observations.
p <- data.frame(peptide_id=c("00123","1E10"),sequence=c("ACDEFGH","ACDEFGH"),
  modifications=c("","1xOxidation [M1]"),modifications_available=TRUE,source_row=2:3,
  annotated_sequence="[K].ACDEFGH.[R]",mapping_status="mapped",match_count=1L,
  matched_accession_count=1L,exclusion_reason="",stringsAsFactors=FALSE)
sites <- data.frame(site_id=c("site1","site2"),record_id=1:2,accession=c("2024-03","MARCH1"),
  bond_left=c(10L,20L),bond_right=c(11L,21L),context_display="ACDE|FGHI",context_complete=TRUE,
  neo_n_count=c(2L,0L),neo_c_count=0L,spanning_count=c(0L,1L),stringsAsFactors=FALSE)
placements <- data.frame(peptide_id=p$peptide_id,record_id=1L,accession="2024-03",start=11L,end=17L,
  protein_length=40L,is_master=TRUE,pd_coordinate_agrees=NA,reported_by_pd=TRUE,left_flank="E",right_flank="I")
evidence <- data.frame(placement_id=c("p1","p2"),peptide_id=p$peptide_id,record_id=1L,accession="2024-03",
  start=11L,end=17L,role="N",bond_left=10L,bond_right=11L,classification="neo_N_candidate",site_id="site1",
  digestion_compatible=FALSE,pre_digest_n_label=NA,is_candidate=TRUE,candidate_reason="non_digestion_compatible")
measurements <- data.frame(peptide_id=rep(p$peptide_id,each=3),measurement_label=rep(c("F1","F2","F3"),2),
  value=c(0,NA,1.25e-12,1.25e12,123.123456789012,-0.125),value_status=c("observed","missing",rep("observed",4)),
  abundance_type="Abundance",source_column=rep(c("=1+1","+SUM(1,2)","@value"),2),
  original_value=c("0","","1.25E-12","1.25E12","123.123456789012","-0.125"),included_in_analysis=TRUE)
merops <- data.frame(site_id="site1",site_accession="2024-03",site_bond_left=10L,site_bond_right=11L,
  site_context="ACDE|FGHI",relation_type=c(rep("reference_site_candidate",3),"motif_compatibility"),
  context_status=c("agreement","mismatch","unknown","agreement"),merops_record_id=c("00123","1E10","03-04","=1+1"),
  protease_id="C01.001",protease_name="Test enzyme",cleavage_type=c("physiological","theoretical","synthetic","non-physiological"),
  reference_raw=c("<%Barr\u00e9 %etal, 2014[20140921A139]%>","<%Wang %etal, Article[20110309A170]%>","unrecognized <%source%>",""),
  reference_ids=c("20140921A139","20110309A170","",""),substrate_accession="2024-03",p1_position="0010",
  substrate_species="Test species",method_code_raw="MS",additional_evidence_code_raw="")
ua <- data.frame(record_id=1L,accession="2024-03",bond_left=10L,bond_right=11L,
  feature_type=c("Signal peptide","Propeptide"),boundary="feature_end",description=c("Predicted signal","Propeptide"),
  evidence_codes=c("ECO:0000255","ECO:0000269"),uniprot_accession="2024-03",source_feature_id=c("f1","f2"),
  feature_key=c("key1","key2"),site_id="site1")
ut <- ua;ut$site_id<-NULL;ut$peptide_id<-p$peptide_id;ut$role<-"N";ut$start<-11L;ut$end<-17L;ut$placement_id<-c("p1","p2")
ur <- data.frame(record_id=1:2,accession=sites$accession,status=c("matched","sequence_mismatch"),
  uniprot_accession=sites$accession,entry_version=c("001","002"),sequence_version="1",
  annotation_date="2024-03-04",sequence_date="2024-03-04",fasta_sequence_sha256="abc",uniprot_sequence_sha256="def")
result <- list(peptides=p,matches=placements,measurements=measurements,cleavage=list(sites=sites,evidence=evidence),
  merops=list(matches=merops),uniprot=list(records=ur,site_annotations=ua,terminal_annotations=ut))
original <- serialize(result,NULL)
tables <- pm_excel_tables(result);summary<-pm_site_summary(result)
stopifnot(identical(original,serialize(result,NULL)),length(tables)==9L,
  nrow(tables[["Cleavage sites"]])==2L,nrow(tables[["Peptides"]])==2L,
  nrow(tables[["Cleavage evidence"]])==2L,nrow(tables[["Abundances"]])==6L,
  identical(tables[["Cleavage evidence"]]$Modifications,p$modifications),
  identical(summary$neo_n_count,sites$neo_n_count),summary$merops_reference_count[1]==3L,
  summary$merops_motif_count[1]==1L,summary$merops_mismatch_count[1]==1L,
  summary$merops_context_unknown_count[1]==1L,summary$merops_reference_count[2]==0L,
  grepl("context mismatch",summary$merops_summary[1]),grepl("sequence differs",summary$uniprot_summary[2]),
  tables[["MEROPS annotations"]]$Reference[1]=="Barr\u00e9 et al., 2014",
  tables[["MEROPS annotations"]]$Reference[2]=="Wang et al., Article",
  tables[["MEROPS annotations"]]$Reference[3]==merops$reference_raw[3])

td<-tempfile("excel_");dir.create(td);path<-file.path(td,"results.xlsx")
metadata<-pm_write_excel(result,path)
stopifnot(identical(metadata$sheets,names(tables)),identical(metadata$rows,vapply(tables,nrow,integer(1))))
read_sheet<-function(sheet) openxlsx::read.xlsx(path,sheet=sheet,check.names=FALSE,sep.names=" ",detectDates=FALSE,skipEmptyRows=FALSE)
xp<-read_sheet("Peptides");xq<-read_sheet("Abundances");xm<-read_sheet("MEROPS annotations")
stopifnot(identical(xp$`Peptide ID`,p$peptide_id),is.numeric(xp$`Input row`),
  identical(xm$`MEROPS record ID`,merops$merops_record_id),
  identical(xm$`Reference P1 position text`,merops$p1_position),
  identical(xm$Reference[1:3],tables[["MEROPS annotations"]]$Reference[1:3]),
  identical(xq$`Source column`,measurements$source_column),is.numeric(xq$Abundance),
  is.na(xq$Abundance[2]),xq$Abundance[1]==0,xq$Abundance[3]>0,
  isTRUE(all.equal(xq$Abundance,measurements$value,tolerance=1e-14)),
  identical(read_sheet("Reference checks")$`Annotation date text`,ur$annotation_date),
  identical(read_sheet("Cleavage sites")$`FASTA record`,c("1","2")))
# Inspect the serialized cell types and formulas, independently of the reader.
unzipped<-file.path(td,"xml");dir.create(unzipped);utils::unzip(path,exdir=unzipped)
worksheets<-list.files(file.path(unzipped,"xl","worksheets"),pattern="\\.xml$",full.names=TRUE)
xml<-lapply(worksheets,function(f) paste(readLines(f,warn=FALSE,encoding="UTF-8"),collapse=""))
stopifnot(!any(vapply(xml,function(x) grepl("<f(?: |>)(?!ont)",x,perl=TRUE),logical(1))))
qtyxml<-xml[[match("Abundances",names(tables))]]
stopifnot(grepl('<c r="A2"[^>]*t="s"',qtyxml),grepl('<c r="C2"[^>]*t="n"[^>]*><v>0</v>',qtyxml),
  !grepl('<c r="C3"[^>]*><v>0</v>',qtyxml))
styles<-paste(readLines(file.path(unzipped,"xl","styles.xml"),warn=FALSE),collapse="")
stopifnot(grepl('numFmtId="49"|formatCode="@"',styles),grepl('state="frozen"',xml[[1]]),grepl('<autoFilter',xml[[1]]))
stopifnot(inherits(try(pm_write_excel(result,path),silent=TRUE),"try-error"))

# New analysis sheets must preserve decimal effects/tiny P values in both the
# stored numeric cell and its display format. Integer counts alone use "0".
analysis <- result
analysis$analysis_tables <- list(pair_measurement_results=data.frame(
  measurement_pair_id=c("001","1E10","missing"),log2_effect=c(0.3577,-0.00003,NA_real_),
  p_value=c(1.25e-12,0.041,NA_real_),adjusted_p_value=c(2.5e-12,0.1,NA_real_),
  n_control=c(3L,4L,NA_integer_),stringsAsFactors=FALSE),
  pair_excluded=data.frame(reason=character(),count=integer()))
analysis_path<-file.path(td,"analysis.xlsx")
analysis_meta<-pm_write_excel(analysis,analysis_path)
analysis_read<-openxlsx::read.xlsx(analysis_path,sheet="Pair tests",check.names=FALSE,skipEmptyRows=FALSE)
expected_analysis<-analysis$analysis_tables$pair_measurement_results
stopifnot(identical(analysis_read$measurement_pair_id,expected_analysis$measurement_pair_id),
  isTRUE(all.equal(analysis_read$log2_effect,expected_analysis$log2_effect,tolerance=1e-14)),
  isTRUE(all.equal(analysis_read$p_value,expected_analysis$p_value,tolerance=1e-14)),
  is.na(analysis_read$p_value[3]),analysis_read$p_value[1]>0)
analysis_xml_dir<-file.path(td,"analysis_xml");dir.create(analysis_xml_dir)
utils::unzip(analysis_path,exdir=analysis_xml_dir)
sheet_xml<-function(name) paste(readLines(file.path(analysis_xml_dir,"xl","worksheets",
  paste0("sheet",match(name,analysis_meta$sheets),".xml")),warn=FALSE),collapse="")
analysis_xml<-sheet_xml("Pair tests")
analysis_styles<-paste(readLines(file.path(analysis_xml_dir,"xl","styles.xml"),warn=FALSE),collapse="")
cell_xfs<-regmatches(analysis_styles,regexpr("<cellXfs[^>]*>.*?</cellXfs>",analysis_styles,perl=TRUE))
xfs<-regmatches(cell_xfs,gregexpr("<xf\\b[^>]*>",cell_xfs,perl=TRUE))[[1]]
cell_xml<-function(address) {
  expression<-paste0('<c\\b[^>]*r="',address,'"[^>]*(?:/>|>.*?</c>)')
  regmatches(analysis_xml,regexpr(expression,analysis_xml,perl=TRUE))
}
attribute<-function(text,name) {
  hit<-regmatches(text,regexec(paste0('\\b',name,'="([^"]*)"'),text,perl=TRUE))[[1]]
  if(length(hit))hit[2] else NA_character_
}
cell_format<-function(address) {
  style<-attribute(cell_xml(address),"s");if(is.na(style))style<-"0"
  format_id<-attribute(xfs[as.integer(style)+1L],"numFmtId")
  if(format_id=="0")return("General")
  if(format_id=="1")return("0")
  custom<-regmatches(analysis_styles,gregexpr("<numFmt\\b[^>]*/>",analysis_styles,perl=TRUE))[[1]]
  matches<-custom[vapply(custom,function(x)identical(attribute(x,"numFmtId"),format_id),logical(1))]
  if(length(matches)!=1L)stop("Unexpected numeric cell format.")
  attribute(matches,"formatCode")
}
cell_value<-function(address) {
  hit<-regmatches(cell_xml(address),regexec("<v>([^<]*)</v>",cell_xml(address)))[[1]]
  if(length(hit))as.numeric(hit[2]) else NA_real_
}
stopifnot(tolower(cell_format("B2"))=="general",tolower(cell_format("C2"))=="general",
  tolower(cell_format("D2"))=="general",cell_format("E2")=="0",
  identical(attribute(cell_xml("B2"),"t"),"n"),identical(attribute(cell_xml("C2"),"t"),"n"),
  isTRUE(all.equal(cell_value("B2"),0.3577,tolerance=1e-14)),
  isTRUE(all.equal(cell_value("C2"),1.25e-12,tolerance=1e-14)),
  cell_value("C2")>0,is.na(cell_value("B4")),is.na(cell_value("C4")),
  grepl('<c r="A1"',sheet_xml("Pair exclusions")),!grepl('<c r="A2"',sheet_xml("Pair exclusions")))

# Optional layers absent are different from a supplied layer with zero hits.
minimal<-result;minimal$merops<-NULL;minimal$uniprot<-NULL
ms<-pm_site_summary(minimal)
stopifnot(all(ms$merops_summary=="Not requested"),all(is.na(ms$merops_reference_count)),
  all(ms$uniprot_status=="not_requested"),length(pm_excel_tables(minimal))==5L)
empty_hits<-result;empty_hits$merops$matches<-merops[FALSE,];empty_hits$uniprot$site_annotations<-ua[FALSE,]
stopifnot(all(pm_site_summary(empty_hits)$merops_reference_count==0L),
  pm_site_summary(empty_hits)$uniprot_summary[1]=="No matching boundary annotation")
pm_write_excel(empty_hits,file.path(td,"no_hits.xlsx"))
empty<-result
empty$cleavage$sites<-sites[FALSE,];empty$cleavage$evidence<-evidence[FALSE,]
empty$peptides<-p[FALSE,];empty$matches<-placements[FALSE,];empty$measurements<-measurements[FALSE,]
empty$merops$matches<-merops[FALSE,];empty$uniprot$site_annotations<-ua[FALSE,]
empty$uniprot$terminal_annotations<-ut[FALSE,];empty$uniprot$records<-ur[FALSE,]
stopifnot(all(vapply(pm_excel_tables(empty),nrow,integer(1))==0L))
pm_write_excel(empty,file.path(td,"empty.xlsx"))
no_cleavage<-minimal;no_cleavage$cleavage<-NULL
stopifnot(nrow(pm_site_summary(no_cleavage))==0L,length(pm_excel_tables(no_cleavage))==3L)

# Invalid lookup keys and unsupported Excel dimensions fail explicitly.
duplicate<-result;duplicate$peptides$peptide_id[2]<-duplicate$peptides$peptide_id[1]
stopifnot(inherits(try(pm_excel_tables(duplicate),silent=TRUE),"try-error"))
long<-minimal;long$peptides$modifications[1]<-paste(rep("a",32768),collapse="")
stopifnot(inherits(try(pm_excel_tables(long),silent=TRUE),"try-error"))
large<-list(peptides=data.frame(peptide_id=rep("id",1048576)))
stopifnot(inherits(try(pm_excel_tables(large),silent=TRUE),"try-error"))

# Public export integration: the Excel opt-in adds a workbook without changing CSVs.
source(file.path(root,"R","mapping.R"),encoding="UTF-8")
source(file.path(root,"R","cleavage.R"),encoding="UTF-8")
source(file.path(root,"R","report.R"),encoding="UTF-8")
fasta_path<-file.path(td,"input.fasta")
writeLines(c(">sp|P1|TEST GN=TEST1","MAAAAPEPTIDER"),fasta_path)
pd<-data.frame("Annotated Sequence"=c("AAA","PEPTIDE"),"Protein Accessions"="P1",
  "Modifications"=c("","1xTest [P1]"),"Abundance: F1: Sample"=c(0,NA),check.names=FALSE)
mapped<-pm_map(pd,pm_read_fasta(fasta_path))
mapped<-pm_build_cleavage(mapped,list(enzyme="trypsin",lysine_blocked=FALSE),digest::digest(file=fasta_path,algo="sha256"))
with_excel<-file.path(td,"with_excel");without_excel<-file.path(td,"without_excel")
one<-pm_export(mapped,with_excel,file.path(root,"inst","templates","report_template.html"),list(fasta=fasta_path),excel=TRUE)
two<-pm_export(mapped,without_excel,file.path(root,"inst","templates","report_template.html"),list(fasta=fasta_path),excel=FALSE)
stopifnot(isTRUE(one$excel_export),identical(two$excel_export,FALSE),
  file.exists(file.path(with_excel,"proteomaster_results.xlsx")),
  !file.exists(file.path(without_excel,"proteomaster_results.xlsx")))
csv<-list.files(with_excel,pattern="\\.csv$")
stopifnot(identical(csv,list.files(without_excel,pattern="\\.csv$")),
  all(vapply(csv,function(name) identical(readBin(file.path(with_excel,name),"raw",file.info(file.path(with_excel,name))$size),
    readBin(file.path(without_excel,name),"raw",file.info(file.path(without_excel,name))$size)),logical(1))))
cat("Excel export tests passed.\n")

# Reader-facing formats suppress numerical display noise without changing cells.
reader<-list(Results=data.frame(`Neo log2FC`=c(3.331e-16,-3.331e-16,0.3577),
  `P-value`=c(1.25e-12,0.04,NA_real_),`BH-adjusted P`=c(2.5e-12,0.08,NA_real_),
  `WT measured (of 6)`=c(3L,6L,0L),`Neo form ID`=c("001","1E10","=1+1"),check.names=FALSE),
  `__Results`=data.frame(`Export key`='["comparison","N","001","002"]',check.names=FALSE))
reader_path<-file.path(td,"reader.xlsx");pm_write_excel_tables(reader,reader_path)
rr<-readxl::read_xlsx(reader_path,sheet="Results")
stopifnot(isTRUE(all.equal(rr$`Neo log2FC`,reader$Results$`Neo log2FC`,tolerance=1e-14)),
  identical(rr$`Neo form ID`,reader$Results$`Neo form ID`),rr$`Neo log2FC`[1]==3.331e-16,rr$`Neo log2FC`[2]== -3.331e-16,rr$`P-value`[1]==1.25e-12,
  identical(openxlsx::sheetVisibility(openxlsx::loadWorkbook(reader_path)),c("visible","hidden")))
reader_xml<-file.path(td,"reader_xml");dir.create(reader_xml);utils::unzip(reader_path,exdir=reader_xml)
rs<-paste(readLines(file.path(reader_xml,"xl","styles.xml"),warn=FALSE),collapse="")
stopifnot(grepl('0.00E\\+00',rs),grepl('0.0005',rs),grepl('0.000',rs))
cat("Reader workbook display precision and hidden export catalog passed.\n")

# Ordinary spreadsheet readers must see real dimensions and no dangling OOXML parts.
reader_sheet<-paste(readLines(file.path(reader_xml,"xl","worksheets","sheet1.xml"),warn=FALSE),collapse="")
stopifnot(grepl('<dimension ref="A1:E4"/>',reader_sheet,fixed=TRUE),grepl('&quot;0.000&quot;',rs,fixed=TRUE))
for(relpath in list.files(reader_xml,pattern="\\.rels$",recursive=TRUE,full.names=TRUE)){
  rel<-paste(readLines(relpath,warn=FALSE),collapse="")
  tags<-regmatches(rel,gregexpr('<Relationship\\b[^>]*/>',rel,perl=TRUE))[[1]]
  base<-if(basename(relpath)==".rels")reader_xml else dirname(dirname(relpath))
  for(tag in tags){
    if(identical(attribute(tag,"TargetMode"),"External"))next
    stopifnot(file.exists(file.path(base,attribute(tag,"Target"))))
  }
}
cat("Workbook dimensions and OOXML relationship reachability passed.\n")
