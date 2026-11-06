# No contrast must keep the workbook's canonical rows available for HTML export.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash="/")))
api <- new.env(parent=globalenv())
for(f in list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))sys.source(f,api)
args <- commandArgs(trailingOnly=TRUE)
td <- if(length(args))args[1] else tempfile("no_comparison_")
cfg <- api$pm_wizard_demo(td);cfg$excel <- TRUE
cfg$sample_design$comparisons <- NULL
cfg$sample_design$samples$Condition <- ""
api$pm_run(cfg,td,file.path(root,"inst","templates","report_template.html"))
out <- file.path(td,"results")
html <- paste(readLines(file.path(out,"report.html"),encoding="UTF-8",warn=FALSE),collapse="\n")
json <- strsplit(strsplit(html,'<script type="application/json" id="data">',fixed=TRUE)[[1]][2],"</script>",fixed=TRUE)[[1]][1]
payload <- api$pm_report_decode_json(json)
prep <- api$pm_report_pair_summaries(payload)
stopifnot(nrow(prep)>0L,all(is.na(prep$comparison)),all(is.na(prep$log2_fc)))
x <- as.data.frame(readxl::read_excel(file.path(out,"atelente_results.xlsx"),sheet="Results"))
stopifnot(nrow(x)>0L,all(x$Status=="No comparison"),all(is.na(x[["Ratio log2FC"]])),
 all(is.na(x[["P-value"]])),all(is.na(x[["BH-adjusted P"]])))
tabs <- api$pm_export_tables(payload);ov <- attr(tabs,"overview")
stopifnot(nrow(ov)==nrow(x),all(ov$Status=="No comparison"))
keys <- lapply(ov[["Export key"]],jsonlite::fromJSON)
stopifnot(all(vapply(keys,function(k)is.na(k[1]),logical(1))))
stopifnot(!any(grepl("^NA measured|^ measured",names(x))))
# Two named conditions without a selected contrast must still not invent one.
cfg$output_dir <- "named";cfg$sample_design$samples$Condition <- rep(c("A","B"),each=3)
api$pm_run(cfg,td,file.path(root,"inst","templates","report_template.html"))
y <- as.data.frame(readxl::read_excel(file.path(td,"named","atelente_results.xlsx"),sheet="Results"))
stopifnot(identical(x,y))
node <- Sys.which("node");stopifnot(nzchar(node))
log <- system2(node,shQuote(c(file.path(root,"tests/test_report_dom_smoke.cjs"),
  file.path(out,"report.html"))),stdout=TRUE,stderr=TRUE)
stopifnot(is.null(attr(log,"status")) || attr(log,"status")==0L)
export_line <- log[startsWith(log,'{"filtered_export":')]
stopifnot(length(export_line)==1L)
export <- jsonlite::fromJSON(export_line)
current <- as.data.frame(readxl::read_excel(file.path(export$filtered_export,"filtered.xlsx"),sheet="Results"))
stopifnot(export$rows==nrow(x),setequal(names(current),names(x)))
order_x <- order(x$`Neo form ID`,x$`Spanning form ID`)
order_current <- order(current$`Neo form ID`,current$`Spanning form ID`)
stopifnot(isTRUE(all.equal(x[order_x,,drop=FALSE],current[order_current,names(x),drop=FALSE],check.attributes=FALSE)))
cat("No-comparison workbook rows, null export keys, absent effects and retained observations passed.\n")
cat("REPORT=",normalizePath(file.path(out,"report.html"),winslash="/"),"\n",sep="")
