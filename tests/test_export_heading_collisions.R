# Non-ASCII literals below must be read as UTF-8 even when Rscript starts in LC_ALL=C.
if(!isTRUE(l10n_info()[["UTF-8"]]))for(loc in c("C.UTF-8","C.utf8","en_US.UTF-8"))if(nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE",loc))))break
args<-commandArgs();root<-dirname(dirname(normalizePath(sub("^--file=","",grep("^--file=",args,value=TRUE)[1]),winslash="/")))
source(file.path(root,"R/report.R"),encoding="UTF-8")
source(file.path(root,"R/excel_locale.R"));source(file.path(root,"R/excel.R"));source(file.path(root,"R/export_bundle.R"))

# Independent fixture: two control neo quantities 2,8; treatment 32. The
# spanning quantities 2,2,4 give ratio log2FC = 3 - mean(0,2) = 2.
samples<-data.frame(Sample=c("İlaç α & <1>","零 🧪","A quoted","excluded"),Condition=c("Base","Base","Dose","Dose"),Use=c("Yes","Yes","Yes","No"),"Peptide column"=paste0("q",1:4),check.names=FALSE)
sites<-data.frame(site_id="site:001",record_id=1L,accession="00123",bond_left=10L,bond_right=11L,context_display="AAAA | GGGG")
peptides<-data.frame(peptide_id=c("00001","spanning","edge"),sequence=c("GGGG","AAAAGGGG","AAAA"),modifications=c("α & <literal>","",NA),source_row=1:3)
evidence<-data.frame(sites[rep(1L,3),],role=c("N","spanning","N"),peptide_id=peptides$peptide_id,start=c(11L,7L,11L),end=14L,is_candidate=c(TRUE,FALSE,TRUE))
terms<-do.call(rbind,lapply(c("N","C"),function(role){x<-evidence;x$role<-role;x$processing_status<-"internal";x$pre_digest_n_label<-FALSE;x$digestion_compatible<-TRUE;x}))
ratios<-data.frame(sites[rep(1L,3),],role="N",neo_peptide_id="00001",spanning_peptide_id="spanning",source_column=paste0("q",1:3),neo_value=c(2,8,32),spanning_value=c(2,2,4),neo_status="observed",spanning_status="observed",log2_ratio=c(0,2,3))
measurements<-do.call(rbind,lapply(1:3,function(i)data.frame(peptide_id=peptides$peptide_id[i],source_column=paste0("q",1:4),value=list(c(2,8,32,2^100),c(2,2,4,2^-100),c(0,NA,1.23456789012345e-120,pi))[[i]],value_status=if(i==3)c("observed","missing","observed","observed") else rep("observed",4))))
d<-list(peptides=peptides,cleavage=list(sites=sites,evidence=evidence,termini=terms,profile=list(enzyme="trypsin",lysine_blocked=FALSE)),pair_view=list(ratios=ratios,measurements=measurements,samples=samples,comparisons=data.frame(Name="Custom ↓",Control="Base",Treatment="Dose")),column_mapping_audit=list(quantities=data.frame(canonical_column=paste0("q",1:4),source_column=c("Original α","raw <2>","third source","=raw4"),input_scale="linear")))
base<-pm_export_tables(d)
stopifnot(base$Results[["Ratio log2FC"]]==2,
  all(c("Custom ↓ peptide log2FC","Custom ↓ Base measured (of 2)","Custom ↓ Dose measured (of 1)",paste0("Quantity: ",samples$Sample)) %in% names(base$Peptides)),
  base$Peptides[["Quantity: İlaç α & <1>"]][3]==0,is.na(base$Peptides[["Quantity: 零 🧪"]][3]),
  base$Peptides[["Quantity: A quoted"]][3]==1.23456789012345e-120)

# Valid custom contrasts previously overwrote a control count (0 with 1).
custom<-d;custom$pair_view$samples$Condition<-c("B C","C","D","D")
custom$pair_view$comparisons<-data.frame(Name=c("A","A B"),Control=c("B C","C"),Treatment=c("D","D"))
j<-custom$pair_view$measurements$peptide_id=="edge"&custom$pair_view$measurements$source_column=="q2"
custom$pair_view$measurements$value[j]<-2;custom$pair_view$measurements$value_status[j]<-"observed"
ct<-pm_export_tables(custom)
stopifnot(!anyDuplicated(names(ct$Peptides)),length(grep("measured",names(ct$Peptides)))==4L,
  ct$Peptides[["A B C measured (of 1) [comparison A; control B C]"]][3]==0L,
  ct$Peptides[["A B C measured (of 1) [comparison A B; control C]"]][3]==1L,
  all(c("A peptide log2FC","A B peptide log2FC","A D measured (of 1)","A B D measured (of 1)") %in% names(ct$Peptides)))

# Duplicate excluded name plus literal suffix. The generated contextual name
# is itself another sample's ordinary name and must be reserved for that sample.
named<-d;named$pair_view$samples$Sample<-c("A","A [Original α]","A [Original α] [source Original α]","A")
stopifnot(!anyDuplicated(named$pair_view$samples$Sample[named$pair_view$samples$Use=="Yes"]))
nt<-pm_export_tables(named)
quantity_labels<-c("Quantity: A [Original α] [source Original α] (2)","Quantity: A [Original α] [source raw <2>]","Quantity: A [Original α] [source Original α]","Quantity: A [=raw4]")
stopifnot(!anyDuplicated(names(nt$Peptides)),all(quantity_labels %in% names(nt$Peptides)),
  identical(nt$Results,base$Results),identical(nt,pm_export_tables(named)))
for(i in seq_along(quantity_labels))stopifnot(identical(nt$Peptides[[quantity_labels[i]]],base$Peptides[[paste0("Quantity: ",samples$Sample[i])]]))
design<-nt$Settings[nt$Settings$Section=="Samples",]
stopifnot(identical(design$Sample,named$pair_view$samples$Sample),identical(design$Include,samples$Use),identical(design[["Source column"]],d$column_mapping_audit$quantities$source_column))

# A quantity label can also collide with comparison effects and coverage.
cross<-d;cross$pair_view$comparisons$Name<-"Quantity: α"
cross$pair_view$samples$Sample[1:2]<-c("α peptide log2FC","α Base measured (of 2)")
xt<-pm_export_tables(cross)
stopifnot(!anyDuplicated(names(xt$Peptides)),
  xt$Peptides[["Quantity: α peptide log2FC [comparison Quantity: α]"]][1]==3,
  xt$Peptides[["Quantity: α peptide log2FC [source Original α]"]][1]==2,
  xt$Peptides[["Quantity: α Base measured (of 2) [comparison Quantity: α; control Base]"]][1]==2L,
  xt$Peptides[["Quantity: α Base measured (of 2) [source raw <2>]"]][1]==8)

# An independent XLSX reader must see both columns and their original values.
output<-tempfile(fileext=".xlsx");pm_write_excel_tables(xt,output,c("Peptides","Settings"))
saved<-as.data.frame(readxl::read_xlsx(output,sheet="Peptides",.name_repair="minimal"))
stopifnot(identical(names(saved),names(xt$Peptides)),!anyDuplicated(names(saved)),
  saved[["Quantity: α peptide log2FC [comparison Quantity: α]"]][1]==3,
  saved[["Quantity: α peptide log2FC [source Original α]"]][1]==2,
  is.na(saved[["Quantity: α Base measured (of 2) [source raw <2>]"]][3]))
cat("Collision-safe headings: custom contrasts, Unicode/excluded names, literal suffixes, effects/counts versus quantities, unchanged values and XLSX round trip passed.\n")
