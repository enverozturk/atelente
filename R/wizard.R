# Explicit column mapping for the local wizard and configuration-driven runs.
# One source row remains one observation; the original table is returned intact.
pm_wizard_frame <- function(x, fields) {
  if (is.null(x) || !length(x)) return(as.data.frame(setNames(rep(list(character()),length(fields)),fields),check.names=FALSE))
  if (!is.data.frame(x)) {
    x <- jsonlite::fromJSON(jsonlite::toJSON(x,auto_unbox=TRUE,null="null"),simplifyDataFrame=TRUE)
    if (!is.data.frame(x)) x <- as.data.frame(x,stringsAsFactors=FALSE,check.names=FALSE)
  }
  if (!setequal(names(x),fields) || anyDuplicated(names(x))) stop("Expected fields: ",paste(fields,collapse=", "))
  x <- x[,fields,drop=FALSE]
  for(n in fields) { x[[n]] <- as.character(x[[n]]);x[[n]][is.na(x[[n]])] <- "" }
  x
}

# Column headers can contain control characters, e.g. wrapped Excel header cells
# exported with "\r\n". A browser does not return such an option value intact
# when the user clicks it, so these columns use a stable ASCII key in the
# selection controls. The server resolves the key back to the exact original
# header; headers without control characters are used as their own key.
pm_wizard_column_keys <- function(columns) {
  columns <- as.character(columns)
  keys <- columns
  unsafe <- which(grepl("[[:cntrl:]]", columns))
  if (length(unsafe)) {
    tokens <- paste0("atelente-column-", unsafe)
    while (any(tokens %in% columns)) tokens <- paste0("_", tokens)  # never equal a real header
    keys[unsafe] <- tokens
  }
  keys
}
pm_wizard_resolve_columns <- function(values, columns) {
  if (is.null(values)) return(values)
  values <- as.character(values)
  hit <- match(values, pm_wizard_column_keys(columns))
  values[!is.na(hit)] <- as.character(columns)[hit[!is.na(hit)]]
  values
}

pm_prepare_column_input <- function(source, column_mapping, sample_design=NULL,
                                    assignment_mode="representative", pair_statistics=FALSE) {
  if(!is.data.frame(source)||anyDuplicated(names(source))||any(!nzchar(names(source))))
    stop("The peptide table needs unique, non-empty column headers.")
  if(!nrow(source)) stop("The peptide table has no observations.")
  if(!is.list(column_mapping)) stop("column_mapping must be a named list.")
  if(!assignment_mode %in% c("unique","representative")) stop("Unknown assignment_mode.")
  if(!is.logical(pair_statistics)||length(pair_statistics)!=1L||is.na(pair_statistics)) stop("pair_statistics must be TRUE or FALSE.")
  selected <- function(name,required=FALSE) {
    x <- column_mapping[[name]]
    if(is.null(x)||identical(x,"")) {if(required)stop("Select the ",name," column.");return(NULL)}
    if(!is.character(x)||length(x)!=1L||is.na(x)||!x %in% names(source))stop("Unknown ",name," column.")
    x
  }
  fields <- list(sequence=selected("sequence",TRUE),modifications=selected("modifications"),
    representative_protein=selected("representative_protein",assignment_mode=="representative"),
    protein_accessions=selected("protein_accessions"),contaminant=selected("contaminant"))
  exclusive <- unlist(fields[c("sequence","modifications","contaminant")],use.names=FALSE)
  proteins <- unlist(fields[c("representative_protein","protein_accessions")],use.names=FALSE)
  if(anyDuplicated(exclusive)||length(intersect(exclusive,proteins)))stop("Sequence, modification and flag mappings must use distinct columns.")
  quantities <- as.character(unlist(column_mapping$quantity_columns,use.names=FALSE))
  if(anyNA(quantities)||anyDuplicated(quantities)||any(!quantities %in% names(source)))stop("Quantity columns must be existing, unique headers.")
  if(length(intersect(quantities,c(exclusive,proteins))))stop("A metadata column cannot also be a quantity column.")
  scale <- column_mapping$quantity_scale
  if(is.null(scale))scale <- "linear"
  if(!is.character(scale)||length(scale)!=1L||!scale %in% c("linear","log2"))stop("quantity_scale must be linear or log2.")
  out <- data.frame(Sequence=as.character(source[[fields$sequence]]),check.names=FALSE,stringsAsFactors=FALSE)
  target <- c(modifications="Modifications",representative_protein="Master Protein Accessions",
    protein_accessions="Protein Accessions",contaminant="Contaminant")
  for(n in names(target))if(!is.null(fields[[n]]))out[[target[[n]]]] <- as.character(source[[fields[[n]]]])
  canonical <- sprintf("Abundance: Selected_%03d",seq_along(quantities))
  conversions <- integer(length(quantities))
  for(i in seq_along(quantities)) {
    values <- as.character(source[[quantities[i]]])
    if(scale=="log2") {
      parsed <- suppressWarnings(as.numeric(trimws(values)))
      finite <- is.finite(parsed)
      linear <- 2^parsed[finite]
      if(any(!is.finite(linear)|linear<=0))stop("Selected log2 values exceed the representable linear range in ",quantities[i],"; no values were clipped.")
      values[finite] <- sprintf("%.17g",linear);conversions[i] <- sum(finite)
    }
    out[[canonical[i]]] <- values
  }
  sf <- c("Sample","Condition","Use","Peptide column","Experimental unit")
  samples <- pm_wizard_frame(sample_design$samples,sf)
  comparisons <- pm_wizard_frame(sample_design$comparisons,c("Name","Treatment","Control"))
  if(length(quantities)&&!nrow(samples))stop("Map each selected quantity column to a sample.")
  if(nrow(samples)) {
    samples$Use <- tolower(trimws(samples$Use))
    if(any(!samples$Use %in% c("yes","no")))stop("Sample Use must be Yes or No.")
    samples$Use <- ifelse(samples$Use=="yes","Yes","No")
    samples$Sample <- trimws(samples$Sample);samples$Condition <- trimws(samples$Condition)
    samples[["Experimental unit"]] <- trimws(samples[["Experimental unit"]])
    cols <- samples[["Peptide column"]]
    if(anyDuplicated(cols)||!setequal(cols,quantities))stop("Sample rows must map each selected quantity column exactly once.")
    active <- samples$Use=="Yes"
    if(!any(active))stop("Include at least one sample.")
    if(any(!nzchar(samples$Sample[active]))||anyDuplicated(samples$Sample[active]))stop("Included sample names must be non-empty and unique.")
    if(nrow(comparisons)) {
      for(n in names(comparisons))comparisons[[n]] <- trimws(comparisons[[n]])
      if(any(!nzchar(as.matrix(comparisons)))||anyDuplicated(comparisons$Name)||
         any(comparisons$Treatment==comparisons$Control)||
         any(!c(comparisons$Treatment,comparisons$Control) %in% samples$Condition[active])||
         any(!nzchar(samples$Condition[active])))stop("Comparisons need a unique name and two different included condition labels.")
      if(anyDuplicated(comparisons[c("Treatment","Control")]))stop("A comparison cannot be repeated.")
    }
    if(pair_statistics) {
      if(!nrow(comparisons))stop("Pair tests require an explicit comparison.")
      units <- samples[["Experimental unit"]][active]
      if(any(!nzchar(units))||anyDuplicated(units))stop("Pair tests require distinct, explicitly declared independent units.")
      for(condition in unique(c(comparisons$Treatment,comparisons$Control)))
        if(sum(active & samples$Condition==condition)<3L)stop("Pair tests require at least three included independent units per compared condition.")
    }
  } else if(nrow(comparisons)||pair_statistics)stop("Comparisons and pair tests require selected quantities and samples.")
  original_samples <- samples
  samples[["Peptide column"]] <- canonical[match(samples[["Peptide column"]],quantities)]
  units <- samples[,c("Sample","Condition","Use","Experimental unit"),drop=FALSE]
  samples[["Experimental unit"]] <- NULL;samples[["Protein column"]] <- rep("",nrow(samples))
  samples <- samples[,c("Sample","Condition","Use","Peptide column","Protein column"),drop=FALSE]
  design <- list(samples=samples,comparisons=comparisons,experimental_units=units)
  audit <- list(schema="atelente-column-mapping-1",row_count=nrow(source),fields=fields,
    quantities=data.frame(source_column=quantities,canonical_column=canonical,input_scale=rep(scale,length(quantities)),
      converted_log2_values=conversions,stringsAsFactors=FALSE),assignment_mode=assignment_mode,
    source_samples=original_samples,source_comparisons=comparisons,
    note="Original rows are not aggregated. Log2 inputs are converted for the linear mapping core; original cells are preserved separately. No missing values are filled.")
  list(peptides=out,original=source,samples=samples,comparisons=comparisons,experimental_units=units,
    design=design,mapping_audit=audit)
}

# Self-contained fixture for UI verification; never a biological benchmark.
pm_wizard_demo <- function(directory) {
  if(dir.exists(directory)&&length(list.files(directory,all.files=TRUE,no..=TRUE)))stop("Demo directory must be empty.")
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  p <- data.frame("Peptide sequence"=c("ACDEK","CDEK","MACDEK","FGHIR","GHIR","FGHIRK"),
    "Modification detail"=c("","","","","",""),"Representative ID"="TEST1",
    "Reported IDs"="TEST1",check.names=FALSE)
  for(i in 1:6)p[[paste0("Signal ",i)]] <- c(100,50,200,80,40,160)*c(1,1.2,0.9,1.6,1.8,1.5)[i]
  pm_write_csv(p,file.path(directory,"peptides.csv"))
  writeLines(c(">sp|TEST1|SYNTHETIC GN=DEMO","MACDEKFGHIRK"),file.path(directory,"database.fasta"))
  samples <- data.frame(Sample=paste0("Sample_",1:6),Condition=rep(c("Control","Treatment"),each=3),Use="Yes",
    "Peptide column"=paste0("Signal ",1:6),"Experimental unit"=paste0("Independent_",1:6),check.names=FALSE)
  config <- list(peptides="peptides.csv",fasta="database.fasta",output_dir="results",excel=FALSE,annotations=FALSE,
    scope="reported",assignment_mode="representative",pair_statistics=FALSE,
    cleavage=list(enzyme="trypsin",lysine_blocked=FALSE,analysis_scope="specificity_v3"),
    column_mapping=list(sequence="Peptide sequence",modifications="Modification detail",representative_protein="Representative ID",
      protein_accessions="Reported IDs",quantity_columns=paste0("Signal ",1:6),quantity_scale="linear"),
    sample_design=list(samples=samples,comparisons=data.frame(Name="Treatment_vs_Control",Treatment="Treatment",Control="Control")),
    metadata=list(notes="Synthetic column-wizard fixture; no biological inference."))
  jsonlite::write_json(config,file.path(directory,"mapping.json"),auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
  invisible(config)
}

pm_wizard_app <- function(run_root=file.path(getwd(),"atelente_runs"), template=NULL) {
  if(!requireNamespace("shiny",quietly=TRUE))stop("Install shiny to use the local column wizard.")
  # Local FASTA and XLSX inputs commonly exceed Shiny's 5 MB default.
  options(shiny.maxRequestSize=max(getOption("shiny.maxRequestSize",5*1024^2),250*1024^2))
  section <- function(title, ...) shiny::div(class="setup-section", shiny::h3(title), ...)
  more <- function(title, ...) shiny::tags$details(class="setup-more",shiny::tags$summary(title),...)
  nav <- function(back, forward, label) shiny::div(class="step-actions",
    if(!is.null(back))shiny::actionButton(back,"Back",class="btn-default"),
    shiny::actionButton(forward,label,class="btn-primary"))
  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$title("Atelente \u00b7 Analysis setup"),shiny::tags$style(shiny::HTML(
      "body{background:#f4f7f8;color:#203c43;font-family:'Segoe UI',Arial,sans-serif;font-size:15px}.container-fluid{max-width:1120px;margin:auto;padding:32px 28px}.brand{display:flex;align-items:center;gap:28px;margin-bottom:28px}.brand h1{display:flex;align-items:center;gap:16px;font-size:44px;letter-spacing:-2px;font-weight:750;margin:0;color:#164e63}.brand .setup-label{border-left:1px solid #c9d8df;padding-left:28px;color:#637782;font-size:16px}.brand-mark{width:34px;height:39px;display:inline-block;flex:0 0 34px;position:relative;border-left:4px solid #246da0;border-right:4px solid #246da0;transform:skew(-13deg);border-radius:3px}.brand-mark:before,.brand-mark:after{content:'';position:absolute;height:4px;background:#246da0;left:-2px;right:-2px;top:13px;transform:rotate(-22deg)}.brand-mark:after{top:24px;width:16px;left:6px;background:#68b4e5}.nav-tabs{display:flex;border:0;gap:8px;margin-bottom:20px}.nav-tabs>li{flex:1;float:none}.nav-tabs>li>a{border:1px solid #dce5ea;border-radius:8px;padding:13px 16px;color:#526b77;background:#fff}.nav-tabs>li.active>a,.nav-tabs>li.active>a:hover,.nav-tabs>li.active>a:focus{background:#164e63;color:white;border-color:#164e63}.tab-content{background:#fff;padding:28px;border:1px solid #dce5ea;border-radius:12px}.step-heading{font-size:24px;font-weight:650;margin:0 0 8px}.step-intro{color:#617780;margin-bottom:24px}.setup-section{margin-bottom:24px}.setup-section h3{font-size:17px;font-weight:650;margin:0 0 16px}.setup-more{border-top:1px solid #e3eaee;padding:14px 0}.setup-more summary::after{content:'+';float:right;font-size:18px;line-height:18px}.setup-more[open]>summary::after{content:'-'} .setup-more summary{cursor:pointer;font-weight:600;color:#466470;margin-bottom:12px}.setup-more:not([open]) summary{margin-bottom:0}.step-actions{display:flex;justify-content:flex-end;gap:10px;border-top:1px solid #e3eaee;padding-top:20px;margin-top:24px}.step-actions .btn:first-child:not(.btn-primary){margin-right:auto}.btn{border-radius:6px;padding:9px 17px}.btn-primary{background:#164e63;border-color:#164e63}.btn-primary:hover,.btn-primary:focus{background:#0e3e51;border-color:#0e3e51}.form-control{border-color:#bcced8;border-radius:5px;min-height:38px;box-shadow:none}.shiny-input-container:not(.shiny-input-container-inline){width:100%;max-width:100%}.control-label{font-weight:600}.sample-table,.table-preview{overflow-x:auto}.sample-row{display:grid;grid-template-columns:minmax(190px,1.5fr) minmax(150px,1fr) minmax(140px,1fr) 60px;min-width:630px;gap:14px;align-items:center;padding:12px 0;border-bottom:1px solid #edf2f4}.sample-row>*{min-width:0}.sample-row .shiny-input-container{width:100%;margin-bottom:0}.sample-row .form-control{width:100%;min-width:0}.source-column{overflow-wrap:anywhere;word-break:normal;color:#59717c;font-size:13px}.units-row{display:grid;grid-template-columns:1fr 1fr;gap:16px;align-items:center;max-width:680px}.units-row .form-group{margin-bottom:8px}.small-note,.help-block{font-size:13px;color:#637983}.summary-box{background:#f3f7fa;padding:18px;border-radius:8px;margin:18px 0}.review-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px}.review-item span{display:block;color:#637983;font-size:12px;margin-bottom:5px}.review-item strong{font-size:15px;overflow-wrap:anywhere}.status{white-space:pre-wrap;overflow-wrap:anywhere;margin-top:16px}.ready-message{color:#28715e}.needs-attention{color:#8c5518}.table-preview table{font-size:13px;white-space:nowrap}.result-actions{display:flex;gap:12px;margin-top:16px}@media(max-width:760px){.container-fluid{padding:16px}.tab-content{padding:18px}.nav-tabs{gap:4px}.nav-tabs>li>a{padding:10px 7px;font-size:12px}.brand h1{font-size:36px}.brand{gap:18px}.brand .setup-label{padding-left:18px}.review-grid{grid-template-columns:1fr 1fr}}"))),
    shiny::div(class="brand",shiny::h1(shiny::span(class="brand-mark",`aria-hidden`="true"),"Atelente"),shiny::span(class="setup-label","Analysis setup")),
    shiny::tabsetPanel(id="step",
      shiny::tabPanel("1 \u00b7 Files",value="files",
        shiny::h2(class="step-heading","Add your files"),shiny::p(class="step-intro","Start with a peptide table and its matching FASTA."),
        shiny::fluidRow(shiny::column(6,shiny::fileInput("peptides","Peptide table",accept=c(".csv",".tsv",".txt",".xlsx"),buttonLabel="Choose file",placeholder="CSV, TSV or Excel"),shiny::uiOutput("sheetUI")),
          shiny::column(6,shiny::fileInput("fasta","Protein sequences",accept=c(".fasta",".fa",".faa",".gz"),buttonLabel="Choose file",placeholder="FASTA or FASTA.gz"))),
        shiny::textOutput("tableInfo"),
        more("Preview table",shiny::div(class="table-preview",shiny::tableOutput("preview"))),
        shiny::p(class="small-note","Files stay in this local R session. Reference downloads send protein accession IDs only."),
        nav(NULL,"toColumns","Next: columns")),
      shiny::tabPanel("2 \u00b7 Columns",value="columns",
        shiny::h2(class="step-heading","Match your columns"),shiny::p(class="step-intro","Choose which columns contain each type of information."),
        shiny::uiOutput("columnsUI"),nav("backFiles","toSamples","Next: samples")),
      shiny::tabPanel("3 \u00b7 Samples",value="samples",
        shiny::h2(class="step-heading","Name your samples"),shiny::p(class="step-intro","One row per selected measurement column. Enter condition names to compare groups."),
        shiny::uiOutput("samplesUI"),shiny::br(),
        section("Comparison",shiny::checkboxInput("compare","Compare two conditions",FALSE),
          shiny::conditionalPanel("input.compare",shiny::fluidRow(
            shiny::column(4,shiny::selectInput("control","Control condition",choices=character())),
            shiny::column(4,shiny::selectInput("treatment","Treatment condition",choices=character())),
            shiny::column(4,shiny::textInput("comparisonName","Comparison name","",placeholder="e.g. Treatment vs Control")))),
          shiny::checkboxInput("pairTests","Include statistical tests",FALSE),
          shiny::conditionalPanel("input.pairTests",
            shiny::p(class="small-note","At least 3 independent samples per condition."),
            shiny::p(class="small-note","Tests assume independent samples. Do not enter technical replicate measurements as separate samples."))),
        nav("backColumns","toRun","Next: review & run")),
      shiny::tabPanel("4 \u00b7 Review & run",value="run",
        shiny::h2(class="step-heading","Review & run"),shiny::p(class="step-intro","Confirm how the samples were prepared, then start."),
        shiny::fluidRow(shiny::column(6,section("Sample preparation",shiny::p("Digestion enzyme: Trypsin"),
          shiny::checkboxInput("blocked","Lysines blocked before digestion",FALSE),
          shiny::conditionalPanel("input.blocked",shiny::textInput("nLabel","N-terminal label name before digestion","",placeholder="Exact modification name")))),
          shiny::column(6,section("Outputs",shiny::checkboxInput("excel","Include Excel workbook",TRUE),
            shiny::checkboxInput("annotations","Include MEROPS & UniProt annotations",TRUE),
            shiny::p(class="small-note","References download once, then use the local cache.")))),
        more("Advanced settings",
          shiny::textInput("runRoot","Save runs in",run_root),
          shiny::conditionalPanel("input.annotations",shiny::checkboxInput("ownReferences","Use my own reference files",FALSE),
            shiny::conditionalPanel("input.ownReferences",shiny::fluidRow(shiny::column(6,shiny::textInput("meropsSnapshot","MEROPS snapshot",""),shiny::textInput("meropsProvenance","MEROPS provenance","")),
              shiny::column(6,shiny::textInput("uniprotSnapshot","UniProt snapshot",""),shiny::textInput("uniprotProvenance","UniProt provenance",""))))),
          more("Cache folder",shiny::textInput("cacheDir","Cache folder","",placeholder="Default user cache"))),
        shiny::uiOutput("reviewCards"),
        shiny::div(class="step-actions",shiny::actionButton("backSamples","Back"),shiny::actionButton("run","Run analysis",class="btn-primary")),
        shiny::div(class="status",shiny::textOutput("status")),shiny::uiOutput("resultsUI"))))
  server <- function(input,output,session) {
    result <- shiny::reactiveVal(NULL);status <- shiny::reactiveVal("Ready.")
    output$sheetUI <- shiny::renderUI({shiny::req(input$peptides)
      if(tolower(tools::file_ext(input$peptides$name))!="xlsx")return(NULL)
      if(!requireNamespace("readxl",quietly=TRUE))return(shiny::p("Install readxl to read XLSX files."))
      sheets<-readxl::excel_sheets(input$peptides$datapath)
      shiny::selectInput("sheet","Worksheet",sheets,selected=sheets[1])})
    source <- shiny::reactive({shiny::req(input$peptides)
      ext <- tolower(tools::file_ext(input$peptides$name));path <- input$peptides$datapath
      # Shiny upload temp names do not retain extensions; pm_read_table dispatches by suffix.
      staged <- file.path(tempdir(),paste0("atelente-upload-",session$token,".",ext))
      if(!file.copy(path,staged,overwrite=TRUE))stop("Cannot stage uploaded table.")
      if(ext=="xlsx")shiny::req(input$sheet)
      pm_read_table(staged,if(ext=="xlsx")input$sheet else NULL)})
    output$preview <- shiny::renderTable({head(source(),8)},striped=TRUE,bordered=FALSE,rownames=FALSE)
    output$tableInfo <- shiny::renderText({x<-source();paste(nrow(x),"rows \u00b7",ncol(x),"columns loaded.")})
    output$columnsUI <- shiny::renderUI({x<-source();choices<-names(x)
      # Labels show the original headers; values are browser-safe keys.
      keyed<-setNames(pm_wizard_column_keys(choices),choices)
      optional<-c("Not supplied"="",keyed)
      shiny::tagList(
        section("Required columns",shiny::fluidRow(
          shiny::column(6,shiny::selectInput("sequence","Peptide sequence",c("Choose a column"="",keyed),selected="")),
          shiny::column(6,shiny::selectInput("representative","Master / representative protein",optional,selected="")))),
        section("Sample measurements",shiny::selectizeInput("quantities","Quantity columns",choices=keyed,selected=character(),multiple=TRUE,
          options=list(placeholder="Choose columns to include")),
          shiny::selectInput("scale","Values are",c("Linear intensities"="linear","Log2-transformed intensities"="log2"),selected="linear"),
          shiny::p(class="small-note","Choose individual samples, not group averages. Leave empty for mapping only.")),
        more("Optional columns",shiny::fluidRow(
          shiny::column(6,shiny::selectInput("modifications","Modifications",optional,selected=""),
            shiny::selectInput("accessions","All reported protein accessions",optional,selected="")),
          shiny::column(6,shiny::selectInput("contaminant","Contaminant flag",optional,selected="")))))})
    # Source-column identity, rather than selection position, owns each input.
    # A new upload/worksheet gets new IDs so delayed browser values cannot leak in.
    sampleGeneration <- 0L
    sampleRegistry <- shiny::reactive({
      columns<-names(source());sampleGeneration<<-sampleGeneration+1L
      state<-new.env(parent=emptyenv())
      state$columns<-columns;state$ids<-paste0(sampleGeneration,"_",seq_along(columns))
      # Default sample names come from headers; line breaks become spaces (a text box cannot hold them).
      state$values<-list(sample=gsub("[[:cntrl:]]+"," ",columns),condition=rep("",length(columns)),use=rep(TRUE,length(columns)))
      state})
    sampleControlId <- function(prefix,column) {
      state<-sampleRegistry();paste0(prefix,"_",state$ids[match(column,state$columns)])
    }
    sampleValues <- shiny::reactive({
      state<-sampleRegistry()
      for(prefix in names(state$values))for(i in seq_along(state$columns)) {
        value<-input[[paste0(prefix,"_",state$ids[i])]]
        if(!is.null(value))state$values[[prefix]][i]<-value
      }
      state$values})
    selectedColumns <- function(values) pm_wizard_resolve_columns(values,names(source()))
    output$samplesUI <- shiny::renderUI({cols<-selectedColumns(input$quantities)
      if(!length(cols))return(shiny::p("No quantity columns selected. Mapping and reference annotation can run without quantities."))
      state<-sampleRegistry();indices<-match(cols,state$columns)
      values<-shiny::isolate(sampleValues())
      rows<-lapply(seq_along(cols),function(i){
        shiny::div(class="sample-row",shiny::span(class="source-column",cols[i]),
          shiny::textInput(sampleControlId("sample",cols[i]),NULL,values$sample[indices[i]],width="100%"),
          shiny::textInput(sampleControlId("condition",cols[i]),NULL,values$condition[indices[i]],width="100%"),
          shiny::checkboxInput(sampleControlId("use",cols[i]),NULL,values$use[indices[i]],width="100%"))})
      shiny::div(class="sample-table",shiny::div(class="sample-row",shiny::strong("Source column"),shiny::strong("Sample"),shiny::strong("Condition"),shiny::strong("Include")),rows)})
    sampleInput <- shiny::reactive({cols<-selectedColumns(input$quantities)
      if(!length(cols))return(pm_wizard_frame(NULL,c("Sample","Condition","Use","Peptide column","Experimental unit")))
      values<-sampleValues();indices<-match(cols,sampleRegistry()$columns)
      data.frame(Sample=values$sample[indices],Condition=values$condition[indices],
        Use=ifelse(values$use[indices],"Yes","No"),
        "Peptide column"=cols,"Experimental unit"=values$sample[indices],check.names=FALSE)})
    shiny::observe({s<-sampleInput();choices<-unique(s$Condition[s$Use=="Yes" & nzchar(s$Condition)])
      c0<-shiny::isolate(input$control);t0<-shiny::isolate(input$treatment)
      if(!length(c0)||!c0 %in% choices)c0<-""
      if(!length(t0)||!t0 %in% choices)t0<-""
      opts<-c("Choose condition"="",setNames(choices,choices))
      shiny::updateSelectInput(session,"control",choices=opts,selected=c0);shiny::updateSelectInput(session,"treatment",choices=opts,selected=t0)})
    for(x in list(c("toColumns","columns"),c("toSamples","samples"),c("toRun","run"),c("backFiles","files"),c("backColumns","columns"),c("backSamples","samples")))local({id<-x[1];tab<-x[2];shiny::observeEvent(input[[id]],shiny::updateTabsetPanel(session,"step",selected=tab))})
    automaticName <- shiny::reactiveVal("")
    shiny::observe({
      if(!isTRUE(input$compare))return()
      control<-input$control;treatment<-input$treatment
      if(is.null(control)||is.null(treatment)||!nzchar(control)||!nzchar(treatment)||control==treatment)return()
      proposed<-paste(treatment,"vs",control)
      current<-input$comparisonName
      previous<-shiny::isolate(automaticName())
      if(is.null(current)||!nzchar(trimws(current))||identical(current,previous)) {
        automaticName(proposed)
        if(!identical(current,proposed))shiny::updateTextInput(session,"comparisonName",value=proposed)
      }
    })
    config <- shiny::reactive({
      if(is.null(input$peptides))stop("Choose a peptide table in step 1.")
      if(is.null(input$fasta))stop("Choose a FASTA file in step 1.")
      if(is.null(input$sequence)||!nzchar(input$sequence))stop("Select the peptide sequence column in step 2.")
      clean<-function(x)if(is.null(x)||!nzchar(x))NULL else x
      column<-function(x)selectedColumns(clean(x))
      mapping<-list(sequence=column(input$sequence),modifications=column(input$modifications),representative_protein=column(input$representative),
        protein_accessions=column(input$accessions),contaminant=column(input$contaminant),quantity_columns=as.list(selectedColumns(input$quantities)),quantity_scale=input$scale)
      comparisonName<-input$comparisonName
      if(isTRUE(input$compare)&&(is.null(clean(input$control))||is.null(clean(input$treatment))))
        stop("Choose both control and treatment conditions in step 3.")
      if(isTRUE(input$compare)&&(is.null(comparisonName)||!nzchar(trimws(comparisonName))))comparisonName<-paste(input$treatment,"vs",input$control)
      comparisons<-if(isTRUE(input$compare))data.frame(Name=comparisonName,Treatment=input$treatment,Control=input$control) else pm_wizard_frame(NULL,c("Name","Treatment","Control"))
      design<-list(samples=sampleInput(),comparisons=comparisons)
      pm_prepare_column_input(source(),mapping,design,"representative",isTRUE(input$pairTests))
      x<-list(peptides=paste0("inputs/peptides.",tolower(tools::file_ext(input$peptides$name))),
        fasta=paste0("inputs/database.",if(grepl("\\.gz$",input$fasta$name,ignore.case=TRUE))"fasta.gz" else "fasta"),
        peptide_sheet=if(tolower(tools::file_ext(input$peptides$name))=="xlsx")input$sheet else NULL,
        output_dir="results",scope="reported",excel=isTRUE(input$excel),annotations=isTRUE(input$annotations),
        assignment_mode="representative",pair_statistics=isTRUE(input$pairTests),column_mapping=mapping,sample_design=design,
        cleavage=list(enzyme="trypsin",lysine_blocked=isTRUE(input$blocked),pre_digest_n_label=if(isTRUE(input$blocked))clean(trimws(input$nLabel)) else NULL,analysis_scope="specificity_v3"),
        cache_dir=clean(trimws(input$cacheDir)),
        metadata=list(sample_unit_policy="Each included quantity column is treated as an independent sample; unit IDs reuse sample names. No pairing or technical-repeat aggregation is inferred."))
      if(isTRUE(input$annotations)&&isTRUE(input$ownReferences)){
      supplied <- 0L
      for(provider in c("merops","uniprot")){
        a<-clean(trimws(input[[paste0(provider,"Snapshot")]]));b<-clean(trimws(input[[paste0(provider,"Provenance")]]))
        if(xor(is.null(a),is.null(b)))stop(provider," requires both snapshot and provenance paths.")
        if(!is.null(a)){x[[provider]]<-list(snapshot=a,provenance=b);supplied<-supplied+1L}
      }
      if(!supplied)stop("Supply a snapshot and its provenance, or turn off Use my own reference files.")
      }
      # Omit unset optional arguments rather than serializing them as JSON null.
      # In particular, an absent cache_dir requests the standard user cache;
      # a named NULL is not a path for pm_run's path resolver.
      x[!vapply(x,is.null,logical(1))]})
    output$reviewCards <- shiny::renderUI({
      if(is.null(input$peptides)||is.null(input$fasta))return(shiny::div(class="summary-box needs-attention","Add both files in step 1."))
      if(is.null(input$sequence)||!nzchar(input$sequence))return(shiny::div(class="summary-box needs-attention","Select a peptide sequence column in step 2."))
      tryCatch({x<-config();item<-function(label,value)shiny::div(class="review-item",shiny::span(label),shiny::strong(value))
        shiny::div(class="summary-box",shiny::div(class="review-grid",
          item("Peptide rows",format(nrow(source()),big.mark=",")),
          item("Included samples",sum(x$sample_design$samples$Use=="Yes")),
          item("Comparison",if(nrow(x$sample_design$comparisons))paste(x$sample_design$comparisons$Treatment,"vs",x$sample_design$comparisons$Control) else "None"),
          item("Digestion",if(x$cleavage$lysine_blocked)"Trypsin / lysines blocked" else "Trypsin"),
          item("Annotations",if(x$annotations)"MEROPS + UniProt" else "Off"),
          item("Statistics",if(x$pair_statistics)"Pair tests" else "Descriptive only")))
      },error=function(e)shiny::div(class="summary-box needs-attention",gsub("representative_protein","master / representative protein",conditionMessage(e),fixed=TRUE)))})
    output$status <- shiny::renderText(status())
    shiny::observeEvent(input$run,{
      result(NULL)
      tryCatch({x<-config();base<-path.expand(trimws(input$runRoot));if(!nzchar(base))stop("Choose an output parent directory.")
        dir.create(base,recursive=TRUE,showWarnings=FALSE);base<-normalizePath(base,winslash="/",mustWork=TRUE)
        job<-tempfile(pattern=paste0("run-",format(Sys.time(),"%Y%m%d-%H%M%S"),"-"),tmpdir=base)
        if(!dir.create(job))stop("Could not create a new run directory.");dir.create(file.path(job,"inputs"))
        stopifnot(file.copy(input$peptides$datapath,file.path(job,x$peptides)),file.copy(input$fasta$datapath,file.path(job,x$fasta)))
        jsonlite::write_json(x,file.path(job,"config.json"),auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
        status("Analysis running. First-use references may take several minutes.")
        shiny::withProgress(message="Running Atelente",value=0.1,{pm_run(x,config_base=job,template=template);shiny::incProgress(0.9)})
        name<-paste0("atelente-result-",gsub("[^A-Za-z0-9]","",session$token))
        shiny::addResourcePath(name,file.path(job,"results"));result(list(directory=job,url=paste0(name,"/report.html")))
        status(paste("Completed. Saved in",job))
      },error=function(e){status(paste("Analysis stopped:",conditionMessage(e)));shiny::showNotification(conditionMessage(e),type="error",duration=NULL)})
    })
    output$resultsUI <- shiny::renderUI({r<-result();if(is.null(r))return(NULL)
      shiny::div(class="result-actions",shiny::tags$a(href=r$url,target="_blank",class="btn btn-primary","Open report"),shiny::downloadButton("reportDownload","Download HTML report"))})
    output$reportDownload <- shiny::downloadHandler(filename=function()"atelente_report.html",content=function(file){r<-result();shiny::req(r);stopifnot(file.copy(file.path(r$directory,"results/report.html"),file))})
  }
  shiny::shinyApp(ui,server)
}
