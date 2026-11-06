args <- commandArgs(trailingOnly=FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=","",grep("^--file=",args,value=TRUE)[1])),".."),winslash="/")
source(file.path(root,"R","report.R"),encoding="UTF-8")
input <- list(label="T\u00fcrk\u00e7e \u03b1", value=pi, zero=0, missing=NULL,
              rows=data.frame(id=seq_len(30000),form=rep("PEPTIDE",30000),value=seq_len(30000)/7))
json <- as.character(jsonlite::toJSON(input,dataframe="rows",auto_unbox=TRUE,digits=NA,null="null"))
packed <- pm_report_encode_json(json)
stopifnot(identical(jsonlite::fromJSON(packed)$schema,"atelente-report-envelope-1"),
          nchar(packed,type="bytes")<nchar(json,type="bytes")/3,
          identical(pm_report_decode_json(packed),jsonlite::fromJSON(json)),
          identical(pm_report_decode_json(json),jsonlite::fromJSON(json)))
temp <- tempfile(fileext=".json")
writeLines(packed,temp,useBytes=TRUE)
# Cross-runtime interchange test; Node is test-only, never required by pm_write_report.
node <- Sys.which("node")
if(nzchar(node)) {
  script <- tempfile(fileext=".cjs")
  js <- sprintf("const fs=require('fs'),C=require(%s);const x=C.parseSync(fs.readFileSync(%s,'utf8'));if(x.rows.length!==30000||x.label!==%s)throw Error('R/browser codec mismatch');",
    jsonlite::toJSON(file.path(root,"inst","templates","report_payload.js"),auto_unbox=TRUE),
    jsonlite::toJSON(temp,auto_unbox=TRUE),jsonlite::toJSON(input$label,auto_unbox=TRUE))
  writeLines(js,script,useBytes=TRUE)
  stopifnot(system2(node,shQuote(script))==0L)
  unlink(script)
}
unlink(temp)
cat("R report transport preserves values and UTF-8; R-to-JavaScript interchange passed.\n")
