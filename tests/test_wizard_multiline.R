# Wizard columns whose headers contain line breaks (wrapped Excel header cells,
# LF or CRLF) must be selectable and must reach the analysis with their exact
# original names. In a browser, clicking an option whose value contains "\r\n"
# does not return that value intact, so such columns use ASCII keys in the
# selection controls and are resolved back on the server.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))source(f,encoding="UTF-8")

# 1. Key and resolution helpers.
cols <- c("Plain header","Wrapped\nLF","Wrapped\r\nCRLF","Tab\there","atelente-column-2")
keys <- pm_wizard_column_keys(cols)
stopifnot(identical(keys[c(1,5)],cols[c(1,5)]),           # ordinary headers are their own key
  !any(grepl("[[:cntrl:]]",keys)),                          # no control characters reach the browser
  !anyDuplicated(keys), !any(keys[2:4] %in% cols),          # keys never collide with a real header
  identical(pm_wizard_resolve_columns(keys,cols),cols),     # keys resolve to the exact headers
  identical(pm_wizard_resolve_columns(cols,cols),cols),     # raw header values remain accepted
  is.null(pm_wizard_resolve_columns(NULL,cols)),
  identical(pm_wizard_resolve_columns(character(),cols),character()),
  identical(pm_wizard_resolve_columns("not a column",cols),"not a column"),
  identical(pm_wizard_column_keys(c("A","B")),c("A","B")))

# 2. Complete wizard run from an XLSX with LF and CRLF headers.
if(requireNamespace("shiny",quietly=TRUE)&&requireNamespace("openxlsx",quietly=TRUE)&&requireNamespace("readxl",quietly=TRUE)&&requireNamespace("writexl",quietly=TRUE)) {
  td <- tempfile("wizard_multiline_");dir.create(td)
  writeLines(c(">sp|P1|ONE","MAARACDEFGHIKLMNRSTVWYK",">sp|P2|TWO","MSTVWYKAGHIKLMNR"),file.path(td,"ref.fasta"))
  table_in <- data.frame(c("CDEFGHIK","DEFGHIK","ACDEFGHIK","CDEFGHI","STVWYK","AGHIK"),
    c("P1","P1","P1","P1","P2","P2"),stringsAsFactors=FALSE)
  names(table_in) <- c("Annotated\nSequence","Master Protein\r\nAccessions")
  quantity <- vapply(1:6,function(i)sprintf(if(i%%2)"Abundance:\nF%d: Sample" else "Abundance:\r\nF%d: Sample",i),character(1))
  for(i in 1:6)table_in[[quantity[i]]] <- c(11,12,13,14,15,16)*1000+i
  expected_headers <- names(table_in)
  # Preserve CRLF as an OOXML-escaped carriage return plus line feed. Some
  # openxlsx versions translate literal XML line breaks on Windows, changing
  # the fixture's headers before the wizard has an opportunity to read them.
  writexl::write_xlsx(list(PeptideGroups=table_in),file.path(td,"peptides.xlsx"))
  stopifnot(identical(names(readxl::read_xlsx(file.path(td,"peptides.xlsx"),sheet="PeptideGroups")),expected_headers))
  app <- pm_wizard_app(run_root=file.path(td,"runs"),template=file.path(root,"inst","templates","report_template.html"))
  shiny::testServer(app$serverFuncSource(), {
    session$setInputs(peptides=list(name="peptides.xlsx",datapath=file.path(td,"peptides.xlsx")),
      fasta=list(name="ref.fasta",datapath=file.path(td,"ref.fasta")),sheet="PeptideGroups")
    headers <- names(source())
    stopifnot(identical(headers,expected_headers))
    key <- setNames(pm_wizard_column_keys(headers),headers)
    # Option values sent to the browser carry no control characters; labels keep the original header text.
    html <- output$columnsUI$html
    option_values <- regmatches(html,gregexpr('value="[^"]*"',html))[[1]]
    stopifnot(length(option_values)>0L,!any(grepl("[\r\n\t]",option_values)),
      all(vapply(paste0('value="',key[grepl("[[:cntrl:]]",headers)],'"'),function(v)v %in% option_values,logical(1))))
    # Browser-returned values are the keys.
    session$setInputs(sequence=unname(key["Annotated\nSequence"]),representative=unname(key["Master Protein\r\nAccessions"]),
      modifications="",accessions="",contaminant="",scale="linear",quantities=unname(key[quantity]),
      compare=FALSE,pairTests=FALSE,ownReferences=FALSE,excel=FALSE,annotations=FALSE,blocked=FALSE,
      cacheDir="",meropsSnapshot="",meropsProvenance="",uniprotSnapshot="",uniprotProvenance="")
    # Default sample names replace header line breaks with a space.
    defaults <- sampleInput()$Sample
    stopifnot(identical(defaults,sprintf("Abundance: F%d: Sample",1:6)))
    for(i in 1:6)do.call(session$setInputs,setNames(list(paste0("S",i)),sampleControlId("sample",quantity[i])))
    s <- sampleInput()
    stopifnot(nrow(s)==6L,identical(s[["Peptide column"]],quantity),identical(s$Sample,paste0("S",1:6)))
    session$setInputs(runRoot=file.path(td,"runs"),run=1L)
    stopifnot(startsWith(status(),"Completed. Saved in"))
    saved <- jsonlite::read_json(file.path(result()$directory,"config.json"),simplifyVector=TRUE)$column_mapping
    stopifnot(identical(saved$sequence,"Annotated\nSequence"),identical(saved$representative_protein,"Master Protein\r\nAccessions"),
      identical(unlist(saved$quantity_columns),quantity),
      file.exists(file.path(result()$directory,"results","report.html")))
  })
}
cat("Wizard multiline headers: LF/CRLF columns selectable by key and passed to the analysis with their exact original names.\n")
