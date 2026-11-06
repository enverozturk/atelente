script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))source(f,encoding="UTF-8")
fails <- function(expr,pattern){e<-tryCatch({force(expr);NULL},error=identity);stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e),fixed=TRUE))}
p <- data.frame(Pep=c("[K].ACDEK.[F]","ACDEK","FGHIR"),Mods=c("","Oxidation [M1]",NA),
  Master=c("P1","P1;P2",""),All=c("P1;P2","P1;P2","P2"),Flag=c("False","","False"),
  First=c("0","2",""),Second=c("-1","1","bad"),Unused=c("\u03b1","001","1e8"),check.names=FALSE)
original <- p
m <- list(sequence="Pep",modifications="Mods",representative_protein="Master",protein_accessions="All",
  contaminant="Flag",quantity_columns=c("First","Second"),quantity_scale="log2")
s <- data.frame(Sample=c("S1","S2"),Condition=c("Control","Treatment"),Use="Yes",
  "Peptide column"=c("First","Second"),"Experimental unit"=c("U1","U2"),check.names=FALSE)
d <- list(samples=s,comparisons=data.frame(Name="T-C",Treatment="Treatment",Control="Control"))
x <- pm_prepare_column_input(p,m,d)
stopifnot(identical(p,original),identical(x$original,original),nrow(x$peptides)==3L,
  identical(x$peptides$Modifications,p$Mods),identical(x$peptides[["Master Protein Accessions"]],p$Master),
  identical(x$peptides[["Abundance: Selected_001"]],c("1","4","")),
  identical(x$peptides[["Abundance: Selected_002"]],c("0.5","2","bad")),
  identical(x$design$samples[["Peptide column"]],c("Abundance: Selected_001","Abundance: Selected_002")),
  identical(x$experimental_units[["Experimental unit"]],c("U1","U2")))
mapped <- pm_map(x$peptides,data.frame(record_id=1L,accession="P1",header="P1",sequence="MACDEKFGHIR",length=10L),sequence_column="Sequence")
stopifnot(nrow(mapped$peptides)==3L,length(unique(mapped$peptides$peptide_id))==3L,
  sum(mapped$measurements$value_status=="missing")==1L,sum(mapped$measurements$value_status=="invalid_numeric_value")==1L)
linear <- m;linear$quantity_scale <- "linear"
xl <- pm_prepare_column_input(p,linear,d)
stopifnot(identical(xl$peptides[["Abundance: Selected_001"]],p$First))
no_mod <- m;no_mod$modifications <- NULL
stopifnot(!"Modifications" %in% names(pm_prepare_column_input(p,no_mod,d)$peptides))
no_quant <- m;no_quant$quantity_columns <- NULL
stopifnot(nrow(pm_prepare_column_input(p,no_quant)$design$samples)==0L)
no_cmp <- d;no_cmp$comparisons <- NULL;no_cmp$samples$Condition <- ""
stopifnot(nrow(pm_prepare_column_input(p,m,no_cmp)$comparisons)==0L)
no_master <- m;no_master$representative_protein <- NULL
fails(pm_prepare_column_input(p,no_master,d),"representative_protein")
stopifnot(!"Master Protein Accessions" %in% names(pm_prepare_column_input(p,no_master,d,"unique")$peptides))
bad <- m;bad$quantity_columns <- c("First","First");fails(pm_prepare_column_input(p,bad,d),"unique headers")
bad <- m;bad$modifications <- "Pep";fails(pm_prepare_column_input(p,bad,d),"distinct columns")
bad <- m;bad$quantity_columns <- "Pep";fails(pm_prepare_column_input(p,bad,d),"metadata column")
bad <- d;bad$samples$Sample <- "same";fails(pm_prepare_column_input(p,m,bad),"unique")
bad <- d;bad$samples[["Peptide column"]] <- "First";fails(pm_prepare_column_input(p,m,bad),"exactly once")
bad <- d;bad$comparisons$Treatment <- "Other";fails(pm_prepare_column_input(p,m,bad),"condition labels")
bad <- p;bad$First[1] <- "1024";fails(pm_prepare_column_input(bad,m,d),"representable linear range")
fails(pm_prepare_column_input(p,m,d,pair_statistics=TRUE),"at least three")
bad <- d;bad$samples[["Experimental unit"]] <- "repeated";fails(pm_prepare_column_input(p,m,bad,pair_statistics=TRUE),"distinct")
fails(pm_prepare_column_input(p,no_quant,pair_statistics=TRUE),"require selected quantities")
td <- tempfile("wizard_demo_");cfg <- pm_wizard_demo(td)
# The annotation-enabled fixture uses a syntactically valid accession. The
# public demo's custom TEST1 identifier is intentionally not queried online.
demo_input <- pm_read_table(file.path(td,"peptides.csv"))
demo_input[["Representative ID"]] <- demo_input[["Reported IDs"]] <- "P00001"
pm_write_csv(demo_input,file.path(td,"peptides.csv"))
writeLines(gsub("TEST1","P00001",readLines(file.path(td,"database.fasta")),fixed=TRUE),file.path(td,"database.fasta"))
loaded <- jsonlite::read_json(file.path(td,"mapping.json"),simplifyVector=TRUE)
input <- pm_read_table(file.path(td,"peptides.csv"))
ready <- pm_prepare_column_input(input,loaded$column_mapping,loaded$sample_design,loaded$assignment_mode,TRUE)
stopifnot(nrow(ready$design$samples)==6L,nrow(ready$design$comparisons)==1L,
  nrow(ready$peptides)==nrow(input),identical(ready$original,input))
if(requireNamespace("shiny",quietly=TRUE)) {
  # Real reference providers/readers with deterministic bytes, no network.
  reference_fields <- c("CLE0000001","A01.001","Test","test","Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly",
    "reference","P00001","4","Test","Test","0","NULL","1-20","MS","","synthetic")
  quoted <- paste0("'",reference_fields,"'");quoted[c(15,18,19)] <- reference_fields[c(15,18,19)]
  reference_raw <- charToRaw(paste0(paste(quoted,collapse="\t"),"\t\r\n"))
  fixture_fetch <- function(url,path) {
    if(grepl("ftp.ebi.ac.uk",url,fixed=TRUE)) {
      writeBin(reference_raw,path);headers<-list("last-modified"="Wed, 22 Feb 2023 12:10:00 GMT")
    } else {
      entry<-list(primaryAccession="P00001",sequence=list(value="MACDEKFGHIRK"),features=list())
      jsonlite::write_json(list(results=list(entry)),path,auto_unbox=TRUE)
      headers<-list("x-uniprot-release"="WIZARD_TEST","x-uniprot-release-date"="1 January 2026","x-total-results"="1")
    }
    list(status_code=200L,headers=headers,url=url,retrieved_utc="2026-11-03T00:00:00Z")
  }
  fixture_cache<-file.path(td,"reference_fixture")
  fixture_merops<-pm_reference_resolve("merops",cache_dir=fixture_cache,fetch=fixture_fetch)
  fixture_uniprot<-pm_reference_resolve("uniprot",cache_dir=fixture_cache,accessions="P00001",fetch=fixture_fetch)
  app <- pm_wizard_app(template=file.path(root,"inst","templates","report_template.html"))
  stopifnot(inherits(app,"shiny.appobj"))
  html <- as.character(app$httpHandler(list(REQUEST_METHOD="GET",PATH_INFO="/",QUERY_STRING=""))$content)
  stopifnot(grepl("Peptide table",html,fixed=TRUE),grepl("Run analysis",html,fixed=TRUE))
  stopifnot(!any(vapply(c('id="preset"','id="saveConfig"','id="autoK"',
    'id="selectionK"','id="declareUnits"'),grepl,logical(1),x=html,fixed=TRUE)))
  shiny::testServer(app$serverFuncSource(), {
    session$setInputs(run=1L)
    stopifnot(grepl("Choose a peptide table",status(),fixed=TRUE))
    session$setInputs(peptides=list(name="peptides.csv",datapath=file.path(td,"peptides.csv")),
      fasta=list(name="database.fasta",datapath=file.path(td,"database.fasta")))
    stopifnot(nrow(source())==6L,grepl("Representative ID",output$columnsUI$html,fixed=TRUE))
    session$setInputs(sequence="Peptide sequence",modifications="Modification detail",
      representative="Representative ID",accessions="Reported IDs",contaminant="",
      scale="linear",quantities=paste0("Signal ",1:6),
      compare=TRUE,comparisonName="Treatment_vs_Control",control="Control",treatment="Treatment",
      pairTests=FALSE,ownReferences=FALSE,excel=FALSE,annotations=FALSE,blocked=FALSE,
      cacheDir="",meropsSnapshot="",meropsProvenance="",uniprotSnapshot="",uniprotProvenance="")
    setSample<-function(i,...) {
      values<-list(...)
      names(values)<-vapply(names(values),sampleControlId,character(1),column=paste0("Signal ",i))
      do.call(session$setInputs,values)
    }
    for(i in 1:6)setSample(i,sample=paste0("Sample_",i),condition=if(i<=3)"Control" else "Treatment",use=TRUE)
    stopifnot(!grepl('id="assignment"',output$columnsUI$html,fixed=TRUE),
      !grepl('id="allFasta"',output$columnsUI$html,fixed=TRUE))
    c0<-config()
    session$setInputs(control=NULL)
    fails(config(),"Choose both control and treatment")
    session$setInputs(control="Control")
    session$setInputs(comparisonName="")
    stopifnot(config()$sample_design$comparisons$Name=="Treatment vs Control")
    session$setInputs(comparisonName="My comparison")
    stopifnot(config()$sample_design$comparisons$Name=="My comparison")
    session$setInputs(comparisonName="Treatment_vs_Control")
    # No duplicate unit form: IDs track sample names, including later edits.
    stopifnot(!grepl("unit_1",output$samplesUI$html,fixed=TRUE),
      grepl("Treatment vs Control",output$reviewCards$html,fixed=TRUE))
    setSample(1,sample="Edited sample");session$setInputs(pairTests=TRUE)
    stopifnot(config()$sample_design$samples$Sample[1]=="Edited sample",
      config()$sample_design$samples[["Experimental unit"]][1]=="Edited sample")
    session$setInputs(pairTests=FALSE)
    stopifnot(config()$sample_design$samples$Sample[1]=="Edited sample")
    setSample(1,sample="Sample_1")
    # Reordering/removing/re-adding selections never reassigns metadata by position.
    # No input values are re-sent to mirror the newly selected order.
    stable_id<-sampleControlId("sample","Signal 1")
    setSample(6,use=FALSE)
    before_reorder<-sampleInput()
    session$setInputs(quantities=paste0("Signal ",6:1))
    reordered<-sampleInput()
    stopifnot(identical(reordered,before_reorder[6:1,,drop=FALSE]) ||
      identical(unname(as.matrix(reordered)),unname(as.matrix(before_reorder[6:1,,drop=FALSE]))),
      identical(sampleControlId("sample","Signal 1"),stable_id),
      grepl('value="Sample_1"',output$samplesUI$html,fixed=TRUE))
    session$setInputs(quantities=paste0("Signal ",2:5))
    stopifnot(identical(sampleInput()$Sample,paste0("Sample_",2:5)))
    # Model unbound browser inputs becoming NULL: cached edits must survive.
    setSample(1,sample=NULL,condition=NULL,use=NULL)
    setSample(6,sample=NULL,condition=NULL,use=NULL)
    session$setInputs(quantities=character())
    stopifnot(nrow(sampleInput())==0L)
    session$setInputs(quantities=paste0("Signal ",1:6))
    stopifnot(identical(sampleInput(),before_reorder),
      grepl('value="Sample_1"',output$samplesUI$html,fixed=TRUE),
      identical(sampleInput()$Use[6],"No"))
    setSample(6,use=TRUE)
    stopifnot(identical(c0$cleavage$analysis_scope,"specificity_v3"),!c0$annotations,
      nrow(c0$sample_design$samples)==6L,identical(c0$column_mapping$sequence,"Peptide sequence"),
      identical(c0$sample_design$comparisons$Treatment,"Treatment"),
      !any(c("cache_dir","peptide_sheet","representative_selection_k") %in% names(c0)),
      !any(vapply(c0,is.null,logical(1))),identical(c0$pair_statistics,FALSE))

    # Trigger the real Run observer: configuration-only tests missed cache_dir=NULL.
    # Blank cache and explicit opt-out must complete a real pm_run/export.
    session$setInputs(excel=TRUE,runRoot=file.path(td,"wizard_runs"),run=3L)
    stopifnot(startsWith(status(),"Completed. Saved in"),!is.null(result()),
      file.exists(file.path(result()$directory,"results","report.html")),
      grepl("Open report",output$resultsUI$html,fixed=TRUE))
    disabled_run<-result()$directory
    book<-file.path(disabled_run,"results","atelente_results.xlsx")
    qp<-openxlsx::read.xlsx(book,sheet="Peptides",check.names=FALSE,sep.names=" ")
    sd<-openxlsx::read.xlsx(book,sheet="Settings",check.names=FALSE,sep.names=" ")
    stopifnot(all(paste0("Quantity: Sample_",1:6) %in% names(qp)),
      identical(sd$`Source column`[sd$Section=="Samples"],paste0("Signal ",1:6)))
    stopifnot(file.exists(file.path(disabled_run,"results","atelente_audit.zip")),
      identical(openxlsx::getSheetNames(file.path(disabled_run,"results","atelente_results.xlsx")),c("Results","Peptides","Settings")))
    saved_config<-jsonlite::read_json(file.path(disabled_run,"config.json"),simplifyVector=TRUE)
    disabled_manifest<-jsonlite::read_json(file.path(disabled_run,"results","run_manifest.json"),simplifyVector=TRUE)
    stopifnot(!"cache_dir" %in% names(saved_config),identical(saved_config$annotations,FALSE),
      identical(disabled_manifest$reference_access$enabled,FALSE))
    # A failed new attempt must not leave a previous successful report presented.
    session$setInputs(representative="",run=4L)
    stopifnot(startsWith(status(),"Analysis stopped:"),is.null(result()),is.null(output$resultsUI))
    session$setInputs(representative="Representative ID")

    # Also exercise enabled annotations with supplied validated snapshots while
    # the cache field stays blank. Any accidental HTTP access fails the test.
    local({
      saved_fetch<-pm_reference_fetch
      on.exit(assign("pm_reference_fetch",saved_fetch,envir=.GlobalEnv))
      assign("pm_reference_fetch",function(...)stop("Unexpected network access in wizard test"),envir=.GlobalEnv)
      session$setInputs(annotations=TRUE,ownReferences=TRUE,meropsSnapshot=fixture_merops$snapshot,meropsProvenance=fixture_merops$provenance,
        uniprotSnapshot=fixture_uniprot$snapshot,uniprotProvenance=fixture_uniprot$provenance,run=2L)
      stopifnot(startsWith(status(),"Completed. Saved in"),!is.null(result()),
        !identical(result()$directory,disabled_run),
        file.exists(file.path(result()$directory,"results","report.html")))
      enabled_manifest<-jsonlite::read_json(file.path(result()$directory,"results","run_manifest.json"),simplifyVector=TRUE)
      stopifnot(isTRUE(enabled_manifest$reference_access$enabled),
        identical(enabled_manifest$reference_access$providers$merops$access,"supplied snapshot"),
        identical(enabled_manifest$reference_access$providers$uniprot$access,"supplied snapshot"),
        identical(enabled_manifest$reference_access$providers$uniprot$release,"WIZARD_TEST"))
    })
    session$setInputs(pairTests=TRUE)
    stopifnot(config()$pair_statistics)
    setSample(1,sample="Sample_2")
    fails(config(),"unique")
    setSample(1,sample="Sample_1",use=FALSE)
    fails(config(),"at least three")
    setSample(1,use=TRUE)
    # Hidden reference paths must not be consumed after switching back to automatic.
    session$setInputs(ownReferences=FALSE)
    stopifnot(!any(c("merops","uniprot") %in% names(config())))
    session$setInputs(ownReferences=TRUE,meropsProvenance="")
    fails(config(),"both snapshot and provenance")
    session$setInputs(annotations=FALSE)
    stopifnot(!any(c("merops","uniprot") %in% names(config())))
    # Removed controls/preset cannot silently override the supported Shiny policy.
    session$setInputs(assignment="unique",allFasta=TRUE,autoK=FALSE,selectionK=99,
      preset=list(name="mapping.json",datapath=file.path(td,"mapping.json")))
    stopifnot(identical(config()$assignment_mode,"representative"),
      identical(config()$scope,"reported"),!"representative_selection_k" %in% names(config()))
    session$setInputs(pairTests=FALSE,compare=FALSE)
    stopifnot(nrow(config()$sample_design$comparisons)==0L)
    session$setInputs(run=5L)
    stopifnot(startsWith(status(),"Completed. Saved in"),!is.null(result()))
    no_contrast <- readxl::read_excel(file.path(result()$directory,"results","atelente_results.xlsx"),sheet="Results")
    stopifnot(nrow(no_contrast)>0L,all(no_contrast$Status=="No comparison"),
      all(is.na(no_contrast$`Ratio log2FC`)))
    session$setInputs(representative="")
    fails(config(),"representative_protein")
    # A different upload with the same headers starts fresh and uses new IDs.
    replacement<-file.path(td,"replacement.csv")
    stopifnot(file.copy(file.path(td,"peptides.csv"),replacement))
    session$setInputs(peptides=list(name="replacement.csv",datapath=replacement))
    stopifnot(!identical(sampleControlId("sample","Signal 1"),stable_id),
      identical(sampleInput()$Sample,paste0("Signal ",1:6)),
      all(sampleInput()$Condition==""),all(sampleInput()$Use=="Yes"))
  })
}
cat("Column wizard mapping, provenance, scales, sample validation, UI construction and reactive workflow tests passed.\n")
