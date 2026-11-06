args <- commandArgs()
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value=TRUE)[1]), winslash="/")))
source(file.path(root, "R", "mapping.R"), encoding="UTF-8")
source(file.path(root, "R", "input_adapter.R"), encoding="UTF-8")
fails <- function(expr) stopifnot(inherits(try(force(expr), silent=TRUE), "try-error"))
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))

fixture <- function() {
  peptides <- data.frame(ID=c("001","002","003","004","005"),
    Sequence=c("[K].PEPTIDE.[R]","PEPTIDE","ACDEFG","ACDEFG","KLMNPQ"),
    Links=c("g1","g1;g1;g2","not-present","","g2"),
    Mods=c("","1xOxidation [M1]",NA,"",""),
    pA=c("8","16","0","-2","NA"), pB=c("4","32","oops","Inf",""),
    Notes=c("keep", "all", "unused", "source", "columns"), check.names=FALSE)
  proteins <- data.frame(Entity=c("g1","g2"),Members=c("P1;P1-2","P2"),
    rA=c("2","4"), rB=c("4","8"), Description=c("first", "second"), check.names=FALSE)
  samples <- data.frame(Sample=c("replicate 1","replicate 2","not selected"),
    Condition=c("same condition","same condition",""), Use=c("Yes","Yes","No"),
    "Peptide column"=c("pA","pB","missing column allowed when unused"),
    "Protein column"=c("rA","rB",""), check.names=FALSE)
  spec <- list(observation_level="peptide_observation",link_separator=";",
    peptide=list(id="ID",sequence="Sequence",protein_links="Links",modifications="Mods",
      scale="linear",quantity_type="reported intensity"),
    protein=list(id="Entity",accessions="Members",entity_type="protein_group",
      scale="linear",quantity_type="reported group intensity"))
  attr(peptides,"source_file") <- "original-peptides.tsv"
  list(peptides=peptides,proteins=proteins,spec=spec,samples=samples)
}
import <- function(x) do.call(pm_import_quant_tables,x)
x <- fixture(); before <- serialize(x,NULL); out <- import(x)
stopifnot(identical(before,serialize(x,NULL)), identical(out$source_peptides,x$peptides),
  identical(out$source_proteins,x$proteins),identical(out$source_samples,x$samples),
  identical(out$spec,x$spec), identical(out$observations$observation_id,x$peptides$ID),
  out$observations$sequence[1]=="PEPTIDE",out$observations$original_sequence[1]=="[K].PEPTIDE.[R]",
  nrow(out$observations)==5L,nrow(out$entities)==2L,nrow(out$links)==5L,
  nrow(out$peptide_measurements)==10L,nrow(out$protein_measurements)==4L,
  identical(out$entities$accessions,x$proteins$Members),all(out$entities$entity_type=="protein_group"),
  all(out$samples$condition=="same condition"),nrow(out$samples)==2L,
  !any(out$observations$excluded),all(is.na(out$observations$contaminant_flag)),
  out$observations$modifications[1]=="",is.na(out$observations$modifications[3]),
  all(out$observations$modifications_available),
  "not-present" %in% out$links$entity_id,!"004" %in% out$links$observation_id)
near(out$peptide_measurements$log2_value[c(1,2,6,7)],c(3,4,2,5))
stopifnot(identical(out$peptide_measurements$value_status,
  c("observed","observed","zero","negative","missing","observed","observed",
    "invalid_numeric_value","invalid_numeric_value","missing")))

# Column spelling/order does not change quantitative observations. No PD master is invented.
y <- x
names(y$peptides) <- c("Row identity","Annotated Sequence","Protein IDs","Modifications",
  "Abundance: F1: Sample","Abundance: F2: Sample","Private note")
names(y$proteins) <- c("Group identity","Protein Accessions","Abundance: F1: Sample",
  "Abundance: F2: Sample","Description")
y$spec$peptide[c("id","sequence","protein_links","modifications")] <-
  list("Row identity","Annotated Sequence","Protein IDs","Modifications")
y$spec$protein[c("id","accessions")] <- list("Group identity","Protein Accessions")
y$samples[["Peptide column"]][1:2] <- y$samples[["Protein column"]][1:2] <-
  c("Abundance: F1: Sample","Abundance: F2: Sample")
y$peptides <- y$peptides[,rev(names(y$peptides)),drop=FALSE]
y$proteins <- y$proteins[,rev(names(y$proteins)),drop=FALSE]
renamed <- import(y)
for (key in c("observations","entities","links")) stopifnot(identical(out[[key]],renamed[[key]]))
for (key in c("peptide_measurements","protein_measurements")) {
  cols <- setdiff(names(out[[key]]),"source_column")
  stopifnot(identical(out[[key]][cols],renamed[[key]][cols]))
}
stopifnot(!any(grepl("master|representative",names(renamed$observations),ignore.case=TRUE)))

# Declared scales preserve zeros, signed log values and extreme finite logs without exponentiation.
v <- pm_input_values(c("0","-2","10000","Inf","-Inf","NaN","1,25","",NA,"NA","2.5e1"),"log2")
stopifnot(identical(v$value_status,c(rep("observed",3),rep("invalid_numeric_value",4),rep("missing",3),"observed")))
near(v$log2_value[c(1:3,11)],c(0,-2,10000,25))
near(v$input_value[c(1:3,11)],c(0,-2,10000,25))
z <- x; z$spec$peptide$scale <- "log2"; z$peptides$pA <- c("0","-2","10000","",NA)
lz <- import(z); near(lz$peptide_measurements$log2_value[1:3],c(0,-2,10000))
stopifnot(all(lz$peptide_measurements$value_status[1:3]=="observed"))

# Missing and blank modification fields are different; IDs remain literal and never become sequence IDs.
z <- x; z$spec$peptide$id <- NULL; z$spec$peptide$modifications <- NULL
generated <- import(z)
stopifnot(identical(generated$observations$observation_id,sprintf("obs_%07d",1:5)),
  !any(generated$observations$modifications_available),all(is.na(generated$observations$modifications)))
z <- x; z$peptides$ID <- c(" 001 ","p1","P1","P1-2","01")
stopifnot(identical(import(z)$observations$observation_id,z$peptides$ID))
z$samples$Condition <- NA_character_
stopifnot(all(import(z)$samples$condition==""))

# Link whitespace and separators are explicit, and groups stay one measured entity.
z <- x; z$peptides$Links[1] <- "g1; g2"; zz <- import(z)
stopifnot(identical(zz$links$entity_id[zz$links$observation_id=="001"],c("g1"," g2")))
z$spec$trim_link_whitespace <- TRUE; zz <- import(z)
stopifnot(identical(zz$links$entity_id[zz$links$observation_id=="001"],c("g1","g2")))
z$spec["link_separator"] <- list(NULL); zz <- import(z)
stopifnot(zz$links$entity_id[zz$links$observation_id=="001"]=="g1; g2")
z <- x; z$peptides$Links[1] <- "P1-2"
stopifnot(import(z)$links$entity_id[1]=="P1-2")

# Only explicit flags exclude; any linked contaminant excludes even a mixed-candidate observation.
z <- x; z$peptides$Flag <- c(FALSE,FALSE,NA,TRUE,FALSE); z$proteins$Flag <- c(FALSE,TRUE)
z$spec$peptide$contaminant <- z$spec$protein$contaminant <- "Flag"
z$peptides$Sequence[4] <- "PEP[+16]TIDE"
zz <- import(z)
stopifnot(identical(zz$observations$excluded,c(FALSE,TRUE,FALSE,TRUE,TRUE)),
  zz$observations$exclusion_reason[2]=="linked_contaminant_entity",
  zz$observations$exclusion_reason[4]=="unsupported_sequence_format;source_contaminant_flag",
  is.na(zz$observations$contaminant_flag[3]),!zz$observations$contaminant_flag[1],
  identical(zz$entities$excluded,c(FALSE,TRUE)),nrow(zz$peptide_measurements)==10L)
z$peptides$Flag <- c("false","false","","+","false"); fails(import(z))
z <- x; z$peptides$Sequence[1] <- ""; zz <- import(z)
stopifnot(zz$observations$excluded[1],zz$observations$exclusion_reason[1]=="missing_sequence")

# Bad mappings/grain/duplicate biological identities fail rather than collapsing rows or guessing.
z <- x; z$peptides$ID[2] <- z$peptides$ID[1]; fails(import(z))
z <- x; z$proteins$Entity[2] <- z$proteins$Entity[1]; fails(import(z))
z <- x; z$peptides$ID[1] <- " "; fails(import(z))
z <- x; z$samples$Sample[2] <- z$samples$Sample[1]; fails(import(z))
z <- x; z$samples[["Peptide column"]][2] <- "pA"; fails(import(z))
z <- x; z$samples[["Protein column"]][2] <- "absent"; fails(import(z))
z <- x; z$spec$peptide$sequence <- "missing"; fails(import(z))
z <- x; z$spec$peptide$scale <- "automatic"; fails(import(z))
z <- x; z$spec$protein$quantity_type <- ""; fails(import(z))
z <- x; z$spec$observation_level <- "PSM"; fails(import(z))
z <- x; z$spec$observation_level <- "feature"; fails(import(z))
z <- x; z$spec$link_separator <- NULL; fails(import(z))
z <- x; z$spec$trim_link_whitespace <- "yes"; fails(import(z))
z <- x; z$samples$Use[1] <- "maybe"; fails(import(z))
z <- x; names(z$peptides)[2] <- names(z$peptides)[1]; fails(import(z))
z <- x; z$peptides$nested <- rep(list(list(1)),nrow(z$peptides)); fails(import(z))

# Empty source tables retain explicit schemas and never synthesize observations/measurements.
z <- x; z$peptides <- z$peptides[FALSE,,drop=FALSE]; z$proteins <- z$proteins[FALSE,,drop=FALSE]
empty <- import(z)
for (key in c("observations","entities","links","peptide_measurements","protein_measurements"))
  stopifnot(nrow(empty[[key]])==0L,identical(names(empty[[key]]),names(out[[key]])))
cat("Input adapter tests passed.\n")
