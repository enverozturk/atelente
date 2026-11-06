script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script, mustWork=TRUE)))
source(file.path(root, "R/report.R"), encoding="UTF-8")
samples <- data.frame(Use="Yes", Condition=c(rep("Control",6),rep("Treatment",5)))
comparisons <- data.frame(Name="T-C",Control="Control",Treatment="Treatment")
x <- list(samples=samples,comparisons=comparisons)
auto <- pm_report_selection(x)
stopifnot(identical(auto$selection$policy,"ends_coverage_overlap_intensity_v3"))
stopifnot(identical(auto$selection$k_by_comparison[["T-C"]],3L))
x$selection_k <- 2
configured <- pm_report_selection(x)
stopifnot(identical(configured$selection$k_by_comparison[["T-C"]],2L),is.null(configured$selection_k))
configured$samples <- samples[1:2,]
stopifnot(identical(pm_report_selection(configured)$selection$k_by_comparison[["T-C"]],2L))
for (invalid in list(0,-1,NA_real_,Inf,1.5,"2",c(2,3))) {
  x$selection_k <- invalid
  stopifnot(inherits(try(pm_report_selection(x),silent=TRUE),"try-error"))
}
x$selection_k <- NULL
x$samples <- samples[c(1,7),]
stopifnot(identical(pm_report_selection(x)$selection$k_by_comparison[["T-C"]],1L))
# Included design only; excluded rows do not increase the frozen default.
x$samples <- samples
x$samples$Use[2:6] <- "No"
stopifnot(identical(pm_report_selection(x)$selection$k_by_comparison[["T-C"]],1L))
serialized <- jsonlite::fromJSON(jsonlite::toJSON(auto,auto_unbox=TRUE,dataframe="rows"))
stopifnot(identical(serialized$selection$k_by_comparison[["T-C"]],3L))
cat("Report analysis selection settings: defaults, fixed overrides, serialization and validation passed.\n")
