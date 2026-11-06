args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]), winslash = "/")))
source(file.path(root, "R", "mapping.R"), encoding = "UTF-8")
source(file.path(root, "R", "report.R"), encoding = "UTF-8")
gene_fixture <- data.frame(record_id=1:3,accession=c("P1","P1-2","P2"),header=c("sp|P1|TEST GN=ABC1 PE=1","sp|P1-2|TEST Isoform 2","sp|P2|TEST no gene"),length=c(10,10,10))
gene_sources <- data.frame(Accession=c("P1","P1-2","P1"),"Gene Symbol"=c("ABC1","ABC1","OLD1"),check.names=FALSE)
gene_result <- pm_report_proteins(gene_fixture,gene_sources)
stopifnot(gene_result$gene_symbol[1]=="ABC1; OLD1",gene_result$gene_symbol[2]=="ABC1",gene_result$gene_symbol[3]=="")
stopifnot(pm_report_proteins(gene_fixture,NULL)$gene_symbol[2]=="")
td <- tempfile("proteomaster_tests_"); dir.create(td)
# Sparse annotations can start beyond readxl's default type-guessing window.
stopifnot(requireNamespace("writexl", quietly = TRUE))
sparse <- data.frame(
  "Annotated Sequence" = rep("AAA", 1005),
  "Master Protein Accessions" = rep("P1", 1005),
  "Modifications in Master Proteins" = c(rep(NA_character_, 1004), "P1 1xAcetyl [N-Term]"),
  "Contaminant" = c(TRUE, FALSE, rep(NA, 1003)),
  "Abundance: F1: Sample" = c(0, 1.25e-12, 1.25e12, rep(NA_real_, 1002)),
  "Source ID" = rep("00123", 1005), check.names = FALSE)
xlsx <- file.path(td, "sparse.xlsx")
writexl::write_xlsx(list(PeptideGroups = sparse), xlsx)
imported <- withCallingHandlers(pm_read_table(xlsx), warning = function(w) stop(w))
stopifnot(nrow(imported) == 1005L,
          imported[["Modifications in Master Proteins"]][1005] == "P1 1xAcetyl [N-Term]",
          all(is.na(imported[["Modifications in Master Proteins"]][1:1004])),
          all(imported[["Source ID"]] == "00123"),
          identical(pm_flag(imported$Contaminant[1:3]), c(TRUE, FALSE, NA)))
im <- pm_measurements(imported)
stopifnot(isTRUE(all.equal(im$value[1:3], c(0, 1.25e-12, 1.25e12))),
          all(im$value_status[1:3] == "observed"), all(im$value_status[4:1005] == "missing"))

# PD repeats intervals without repeating the preceding accession.
positions <- pm_positions("P1 [2-4]; [3-5]; P2 [8-10]; [12-14]")
stopifnot(positions$parsed, identical(positions$data$accession, c("P1", "P1", "P2", "P2")),
          identical(positions$data$start, c(2L,3L,8L,12L)))
stopifnot(pm_positions("P1[2-4]; P2 [8-10]")$parsed,
          nrow(pm_positions("P1 [2-4]; P1 [2-4]")$data) == 2L,
          pm_positions(NA_character_)$parsed, pm_positions("")$parsed)
orphan <- pm_positions("[2-4]; P2 [8-10]")
broken <- pm_positions("P1 [2-4]; unreadable; [3-5]; P2 [8-10]; [12-14]")
stopifnot(!orphan$parsed, identical(orphan$data$accession, "P2"), !broken$parsed,
          identical(broken$data$accession, c("P1", "P2", "P2")),
          !pm_positions("P1 [999999999999999-999999999999999]")$parsed)
fa <- file.path(td, "tiny.fasta")
writeLines(c(">sp|P1|ONE test", "MAAAAKPEPTIDER", ">sp|P1-2|TWO test", "KPEPTIDER", ">sp|P2|THREE test", "KPEPTIDER",
             ">sp|C1|CONTAM test", "MCONTAMR", ">sp|P1|DUPLICATE test", "MAAAAKPEPTIDER"), fa)
f <- pm_read_fasta(fa)
stopifnot(nrow(f) == 5L, sum(f$accession == "P1") == 2L)
gz <- file.path(td, "tiny.fasta.gz"); con <- gzfile(gz, "wt"); writeLines(readLines(fa), con); close(con)
stopifnot(identical(f, pm_read_fasta(gz)))
p <- data.frame("Annotated Sequence" = c("[M].AAA.[A]", "[K].PEPTIDE.[R]", "CONTAM", "MISSING", "AAA", "PEPTIDE", "A[+80]AA"),
 "Protein Accessions" = c("P1", "P1; P1-2", "C1", "ABSENT", "P1", "P1; P2", "P1"),
 "Master Protein Accessions" = c("P1", "P1", "C1", "ABSENT", "P1", "P1", "P1"),
 "Positions in Master Proteins" = c("P1 [2-4]", "P1 [6-12]", "C1 [2-7]", "ABSENT [1-7]", "P1 [2-4]", "P1 [7-13]", ""),
 "Modifications" = c("", "1xTest [P1]", "", "", "", "", ""),
 "Contaminant" = c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE),
 "Abundance: F1: Sample" = c(10,NA,30,0,50,60,70),
 check.names = FALSE)
r <- pm_map(p, f)
# AAA occurs twice, including the overlap, in each duplicate P1 record.
stopifnot(r$peptides$match_count[1] == 4L, identical(sort(unique(r$matches$start[r$matches$peptide_id == "pep_0000001"])), c(2L,3L)))
stopifnot(r$peptides$match_count[2] == 3L, all(c("P1","P1-2") %in% r$matches$accession[r$matches$peptide_id == "pep_0000002"]))
stopifnot(r$peptides$mapping_status[3] == "excluded_contaminant", !"pep_0000003" %in% r$matches$peptide_id)
stopifnot(r$peptides$mapping_status[4] == "unmatched", r$candidates$status[r$candidates$peptide_id == "pep_0000004"] == "accession_absent")
stopifnot(r$peptides$mapping_status[7] == "unsupported_sequence_format", r$peptides$modifications[2] == "1xTest [P1]")
stopifnot(is.na(r$measurements$value[2]), r$measurements$value[4] == 0, !r$measurements$included_in_analysis[3])
stopifnot(r$pd_checks$status[r$pd_checks$peptide_id == "pep_0000002"] == "different_position")
ra <- pm_map(p[2,,drop=FALSE], f, scope = "all")
stopifnot(ra$peptides$match_count == 4L, any(!ra$matches$reported_by_pd), ra$matches$is_master[ra$matches$accession == "P2"] == FALSE)
bare <- data.frame("Annotated Sequence" = "AAA", check.names = FALSE)
rb <- pm_map(bare, f, scope = "all")
stopifnot(!pm_report_peptides(rb$peptides,bare)$modifications_available,
          is.na(pm_report_peptides(rb$peptides,bare)$pd_master_accessions))
report_peptides <- pm_report_peptides(r$peptides,p)
stopifnot(all(report_peptides$modifications_available),report_peptides$pd_master_accessions[1]=="P1",
          report_peptides$modifications[1]=="",report_peptides$modifications[2]=="1xTest [P1]")
stopifnot(rb$peptides$contaminant_screening == "not_available", is.na(rb$peptides$modifications), nrow(rb$measurements) == 0)
stopifnot(inherits(try(pm_map(bare, f), silent = TRUE), "try-error"))
# Pre-cleaned inputs need neither contaminant fields nor a contaminant list.
clean <- data.frame("Annotated Sequence" = "AAA", "Master Protein Accessions" = "P1",
  "Positions in Master Proteins" = "P1 [2-4]; [3-5]", check.names = FALSE)
clean_result <- withCallingHandlers(pm_map(clean, f), warning = function(w) stop(w))
stopifnot(nrow(clean_result$excluded) == 0L, nrow(clean_result$matches) == 4L,
          nrow(clean_result$pd_checks) == 2L, all(clean_result$pd_checks$status == "coordinate_match"),
          all(clean_result$matches$pd_coordinate_agrees), clean_result$peptides$pd_position_format_valid,
          is.na(clean_result$peptides$contaminant_flag))
rc <- pm_map(p[6,,drop=FALSE], f, contaminant_ids = "P2")
stopifnot(rc$excluded$exclusion_reason == "shared_contaminant_assignment")
# Sequence overlap detects an external contaminant even with a different accession.
cf <- data.frame(accession = "External", sequence = "MPEPTIDER")
rc2 <- pm_map(p[2,,drop=FALSE], f, contaminant_fasta = cf)
stopifnot(rc2$excluded$exclusion_reason == "contaminant_sequence_overlap")
prot <- data.frame(Accession = "P1", Contaminant = TRUE)
stopifnot(nrow(pm_map(p[1,,drop=FALSE], f, proteins = prot)$excluded) == 1L)
stopifnot(identical(pm_flag(c("YANLI\u015e", "TRUE", "unknown", "")), c(FALSE, TRUE, NA, NA)))
stopifnot(identical(pm_occurrences("AA", "AAAA"), c(1L,2L,3L)), length(pm_occurrences("IL", "LL")) == 0)
stopifnot(pm_parse_sequence("peptide")$sequence == "PEPTIDE", pm_upper("il") == "IL")
stopifnot(!pm_positions("P1 [2-4]; unreadable")$parsed)
terminal <- data.frame("Annotated Sequence" = "KPEPTIDER", "Protein Accessions" = "P1-2", check.names=FALSE)
rt <- pm_map(terminal, f)
stopifnot(rt$matches$start == 1L, rt$matches$end == 9L, rt$matches$left_flank == "-", rt$matches$right_flank == "-")
# Distinct modified observations keep distinct row IDs and measurements.
mods <- p[rep(2,2),,drop=FALSE]; mods$Modifications <- c("1xTest [P1]", "2xTest [P1; E2]")
rm <- pm_map(mods, f)
stopifnot(length(unique(rm$matches$peptide_id)) == 2L, nrow(rm$measurements) == 2L)
# Unknown text must not become a false contaminant flag; grouped-only input is valid.
unknown <- bare; unknown$Contaminant <- "unrecognized"; unknown[['Abundances (Grouped): G1']] <- "n/a"
ru <- pm_map(unknown, f, scope="all")
stopifnot(ru$peptides$contaminant_screening == "not_available", is.na(ru$measurements$value), ru$measurements$abundance_type == "grouped")
bad <- file.path(td, "bad.fasta"); writeLines(c(">empty", ">next", "A"), bad)
stopifnot(inherits(try(pm_read_fasta(bad), silent = TRUE), "try-error"))
# Report supports all-excluded and no-quantity inputs without dropping audit rows.
ex <- pm_map(p[3,,drop=FALSE], f)
out <- file.path(td,"out")
pm_export(ex, out, file.path(root,"inst","templates","report_template.html"), list(fasta = fa))
stopifnot(file.exists(file.path(out,"report.html")), nrow(read.csv(file.path(out,"excluded_peptides.csv"), fileEncoding="UTF-8-BOM")) == 1)
stopifnot(inherits(try(pm_export(ex, out, file.path(root,"inst","templates","report_template.html"), list(fasta=fa)), silent=TRUE),"try-error"))
out2 <- file.path(td,"bare")
rb$peptides$modifications <- "</script><script>alert(1)</script>"
pm_export(rb, out2, file.path(root,"inst","templates","report_template.html"), list(fasta=fa))
html <- paste(readLines(file.path(out2,"report.html"), warn=FALSE), collapse="\n")
stopifnot(!grepl("</script><script>alert(1)</script>",html,fixed=TRUE), grepl("\\u003c",html,fixed=TRUE))
cat("All mapping tests passed.\n")
