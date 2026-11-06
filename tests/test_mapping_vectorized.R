# Column-wise parsers used by pm_map must equal the scalar pm_split() and
# pm_parse_sequence() element by element, including blanks, NA, whitespace,
# duplicate/empty accession pieces, flanked and unsupported sequences.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/")))
source(file.path(root,"R","mapping.R"),encoding="UTF-8")
accessions <- c(NA,""," ","P1","P1;P2"," P1 ; P2 ;P1","P1;;P2",";",";;"," ; ","P1;","Q9-2; Q9","ç;x")
sequences <- c(NA,""," ","PEPTIDEK","peptidek"," PEPTIDEK ","[K].PEPTIDER.[A]","[-].MAAR.[G]","[K].pep.[R]",
  "PEP-TIDE","PEP TIDE","[K]PEPTIDE[R]","[K].PEPTIDE.[R] ","123","[].PEP.[R]","ÇAB")
set.seed(3)
pool <- c(LETTERS,letters,";"," ","[","]",".","-")
random <- vapply(1:400,function(i) paste(sample(pool,sample(0:12,1),TRUE),collapse=""),character(1))
for(x in list(accessions,c(accessions,random))) {
  got <- pm_split_all(x)
  for(i in seq_along(x)) if(!identical(got[[i]],pm_split(x[i])))
    stop("pm_split_all differs for ",encodeString(x[i],quote='"'))
}
for(x in list(sequences,c(sequences,random))) {
  got <- pm_parse_sequence_all(x)
  for(i in seq_along(x)) {
    want <- pm_parse_sequence(x[i])
    have <- list(sequence=got$sequence[i],left=got$left[i],right=got$right[i],status=got$status[i])
    if(!identical(have,want)) stop("pm_parse_sequence_all differs for ",encodeString(x[i],quote='"'))
  }
}
stopifnot(identical(pm_split_all(character()),list()),all(lengths(pm_parse_sequence_all(character()))==0L))
cat("Column-wise accession splitting and sequence parsing equal the scalar parsers.\n")
