script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
source(file.path(root,"R/report.R"),encoding="UTF-8")
td <- tempfile("large-reference-report-");dir.create(td)
set.seed(441)
# Incompressible enough to remain above substring's former implicit 1e6 cap.
value <- paste0(sample(c(letters,LETTERS,0:9),1500000,replace=TRUE),collapse="")
for(enabled in c(FALSE,TRUE)) {
  payload <- list(cleavage=list(),text=value,precision=data.frame(value=c(pi,1e-100,1+2^-52,0,NA_real_)),reference_access=if(enabled)
    list(enabled=TRUE,providers=list(uniprot=list(provider="uniprot",release="test"))) else list(enabled=FALSE))
  path <- file.path(td,paste0(enabled,".html"))
  pm_write_report(payload,file.path(root,"inst","templates","report_template.html"),path)
  html <- paste(readLines(path,encoding="UTF-8",warn=FALSE),collapse="\n")
  stopifnot(nchar(html)>1000000,grepl("</html>",html,fixed=TRUE))
  marker <- '<script type="application/json" id="data">'
  start <- regexpr(marker,html,fixed=TRUE)[1]+nchar(marker)
  tail <- substring(html,start,nchar(html))
  end <- regexpr("</script>",tail,fixed=TRUE)[1]
  decoded <- pm_report_decode_json(substr(tail,1,end-1),simplifyVector=FALSE)
  stopifnot(end>0,identical(decoded$payload_schema,"atelente-report-data-2"))
  # The writer deliberately adds its transport schema. Every original field,
  # including the full large string and doubles, must remain exactly unchanged.
  decoded$payload_schema <- NULL
  stopifnot(identical(decoded,payload),is.null(payload$payload_schema))
}
cat("Large reference-enabled/disabled reports preserve complete payload and closing scripts.\n")
