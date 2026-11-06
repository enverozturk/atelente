# The vectorized terminus classifier must equal the scalar pm_terminal() for
# every element, across profiles, labels, protein ends, initiator Met and the
# proline rule. pm_build_cleavage relies on this equivalence.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
source(file.path(root,"R","mapping.R"),encoding="UTF-8")
source(file.path(root,"R","cleavage.R"),encoding="UTF-8")
set.seed(7)
aa <- strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
seqs <- c("MKRPAKRK",vapply(1:40,function(i) paste(c(if(i%%2) "M",sample(c(aa,"K","R","P"),sample(3:30,1),TRUE)),collapse=""),character(1)))
profiles <- list(list(enzyme="trypsin",lysine_blocked=FALSE),list(enzyme="trypsin",lysine_blocked=TRUE),list(enzyme="unknown"))
for(profile in profiles) {
  sequence <- character(); k <- integer(); role <- character(); label <- logical()
  for(s in seqs) for(kk in 0:nchar(s)) for(r in c("N","C")) for(l in c(NA,TRUE,FALSE)) {
    sequence <- c(sequence,s); k <- c(k,kk); role <- c(role,r); label <- c(label,l)
  }
  v <- pm_terminal_vec(sequence,k,role,profile,label)
  for(i in seq_along(k)) {
    s <- pm_terminal(sequence[i],k[i],role[i],profile,label[i])
    for(field in names(s)) if(!identical(s[[field]],v[[field]][i]))
      stop(sprintf("pm_terminal_vec differs for %s (enzyme=%s, k=%d, role=%s, label=%s): %s vs %s",field,
        profile$enzyme,k[i],role[i],label[i],format(s[[field]]),format(v[[field]][i])))
  }
  stopifnot(identical(names(v),c("digestion_compatible","pre_digest_n_label","processing_status",
    "classification","is_candidate","candidate_reason")), all(lengths(v)==length(k)))
}
empty <- pm_terminal_vec(character(),integer(),character(),profiles[[1]],logical())
stopifnot(all(lengths(empty)==0L))
# Placement specificity: malformed placements are still rejected.
bad <- data.frame(placement_id=c("p1","p1","p2"),role=c("N","N","C"),digestion_compatible=NA,
  processing_status="not_assessed",peptide_id="x",record_id=1L,accession="A",start=1L,end=2L)
stopifnot(inherits(try(pm_placement_specificity(bad),silent=TRUE),"try-error"))
cat("Vectorized terminus classification equals pm_terminal for all profiles; malformed placements rejected.\n")
