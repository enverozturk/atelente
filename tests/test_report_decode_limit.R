script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-dirname(dirname(normalizePath(script,winslash="/")))
source(file.path(root,"R/report.R"),encoding="UTF-8")
expect_limit<-function(x,limit){
  error<-tryCatch(pm_report_decode_json(x,max_uncompressed_bytes=limit),error=identity)
  stopifnot(inherits(error,"error"),grepl("refresh decode limit",conditionMessage(error),fixed=TRUE))
}
legacy<-paste0('{"text":"',paste(rep("abc",1000),collapse=""),'"}')
expect_limit(legacy,1000)
stopifnot(identical(pm_report_decode_json(legacy,max_uncompressed_bytes=4000),jsonlite::fromJSON(legacy)))
compressed<-pm_report_encode_json(legacy,threshold=0)
expect_limit(compressed,1000)
stopifnot(identical(pm_report_decode_json(compressed,max_uncompressed_bytes=4000),jsonlite::fromJSON(legacy)))
# Invalid compressed bytes deliberately prove the size check runs BEFORE decode.
oversized<-'{"schema":"atelente-report-envelope-1","encoding":"deflate-base64","uncompressed_bytes":1535736050,"data":"not base64"}'
expect_limit(oversized,450*1024^2)
for(limit in list(NA_real_,0,-1,"100",numeric()))
  stopifnot(inherits(try(pm_report_decode_json("{}",max_uncompressed_bytes=limit),silent=TRUE),"try-error"))
cat("Plain/compressed decode limits reject before allocation; bounded legacy reads remain compatible.\n")
