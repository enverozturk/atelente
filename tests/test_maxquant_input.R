script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- normalizePath(file.path(dirname(script), ".."), winslash="/", mustWork=TRUE)
for (module in c("mapping", "input_adapter", "sample_ratios", "maxquant_input"))
  source(file.path(root, "R", paste0(module, ".R")), encoding="UTF-8")
fails <- function(expr) stopifnot(inherits(try(force(expr), silent=TRUE), "try-error"))
near <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
fixture <- function() {
  peptides <- data.frame(id=c("0","01","2","3","4","5","6","7"),
    Sequence=c("PEPTIDE","PEPTIDE","ACDEFG","KLMNPQ","PEPTIDE","ACDEFG","KLMNPQ","ACDEFG"),
    "Protein group IDs"=c("0","0;01","2;0","3;0","unknown","","0","0"),
    "Mod. peptide IDs"=c("0;1","2","3","4","5","6","7","8"),
    "Potential contaminant"=c("","","","","","","+",""),
    Reverse=c("","","","","","","","+"),
    "Intensity A"=c("8","16","32","64","8","16","32","64"),
    "Intensity B"=c("0","","32","64","8","16","32","64"),
    Intensity="1000", Note="Preserve source summaries", check.names=FALSE)
  proteins <- data.frame(id=c("0","01","2","3"),"Protein IDs"=c("P1;P1-2","P2","P3","P4"),
    "Potential contaminant"=c("","+","",""),Reverse=c("","","+",""),
    "Only identified by site"=c("","","","+"),
    "Intensity A"=c("2","4","8","16"),"Intensity B"=c("4","8","16","32"),
    "LFQ intensity A"=c("1","2","4","8"),"LFQ intensity B"=c("2","4","8","16"),
    "iBAQ A"=c("0.5","1","2","4"),"iBAQ B"=c("1","2","4","8"),
    Intensity="9999",check.names=FALSE)
  samples <- data.frame(Sample=c("01","02"),Condition="",Use="Yes",
    "Peptide column"=c("Intensity A","Intensity B"),"Protein column"=c("Intensity A","Intensity B"),check.names=FALSE)
  attr(peptides,"source_file") <- "original/peptides.txt"
  list(peptides=peptides,proteins=proteins,samples=samples)
}
run <- function(x,...) do.call(pm_import_maxquant_tables,c(x,list(...)))
x <- fixture();before <- serialize(x,NULL);data <- run(x);result <- pm_sample_ratios(data)
stopifnot(identical(before,serialize(x,NULL)),identical(data$source_peptides,x$peptides),
  identical(data$source_proteins,x$proteins),identical(data$source_samples,x$samples),
  identical(data$observations$observation_id,x$peptides$id),identical(data$entities$entity_id,x$proteins$id),
  all(data$observations$source_grain=="sequence_summary"),data$spec$source_grain=="sequence_summary",
  all(is.na(data$observations$modifications)),!any(data$observations$modifications_available),
  identical(data$entities$accessions,x$proteins[["Protein IDs"]]),all(data$entities$entity_type=="protein_group"),
  nrow(data$peptide_measurements)==16L,nrow(data$protein_measurements)==8L,
  !any(grepl("master|representative",names(data$observations),ignore.case=TRUE)),
  identical(data$samples$source_experiment,c("A","B")))
# Linked contaminants/decoys exclude the whole observation, including clean candidate baselines.
stopifnot(identical(data$observations$excluded,c(FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE)),
  data$observations$exclusion_reason[2]=="linked_contaminant_entity",
  data$observations$exclusion_reason[3]=="linked_decoy_entity",
  data$observations$exclusion_reason[7]=="source_contaminant_flag",
  data$observations$exclusion_reason[8]=="source_decoy_flag",
  all(result$ratios$ratio_status[result$ratios$observation_id %in% c("01","2","6","7")]=="excluded_observation"))
# Site-only groups exclude their denominator; another explicitly linked clean group remains usable.
rows <- result$ratios[result$ratios$observation_id=="3",]
stopifnot(all(rows$ratio_status[rows$entity_id=="3"]=="excluded_protein_entity"),
  all(rows$ratio_status[rows$entity_id=="0"]=="normalized"),
  data$entities$only_identified_by_site_flag[4],data$entities$exclusion_reason[4]=="only_identified_by_site")
allowed <- run(x,exclude_only_identified_by_site=FALSE)
stopifnot(allowed$entities$only_identified_by_site_flag[4],!allowed$entities$excluded[4],
  !allowed$spec$exclude_only_identified_by_site)
stopifnot(all(pm_sample_ratios(allowed)$ratios$ratio_status[pm_sample_ratios(allowed)$ratios$observation_id=="3"]=="normalized"))
stopifnot(all(result$ratios$ratio_status[result$ratios$observation_id=="4"]=="protein_entity_not_found"),
  all(result$ratios$ratio_status[result$ratios$observation_id=="5"]=="missing_protein_link"))
row <- result$ratios[result$ratios$observation_id=="0",]
near(row$normalized_ratio,c(4,NA));stopifnot(identical(row$ratio_status,c("normalized","peptide_zero")))
stopifnot(data$peptide_measurements$value_status[10]=="missing",data$peptide_measurements$input_value[9]==0)

# Alternative denominator families are explicit and tied to the same exact experiment suffix.
for (family in c("LFQ intensity","iBAQ")) {
  z <- x;z$samples[["Protein column"]] <- paste(family,c("A","B"));d <- run(z,protein_quantity=family)
  ratio <- pm_sample_ratios(d)$ratios;ratio <- ratio[ratio$observation_id=="0" & ratio$sample_id=="01",]
  near(ratio$normalized_ratio,if(family=="LFQ intensity") 8 else 16)
  stopifnot(d$spec$protein$quantity_type==paste("MaxQuant protein-group",family),
    d$spec$peptide$scale=="linear",d$spec$protein$scale=="linear")
}
z <- x;z$samples[["Protein column"]] <- rev(z$samples[["Protein column"]]);fails(run(z))
z <- x;z$samples <- z$samples[1,,drop=FALSE];z$samples[["Peptide column"]] <- "Intensity";fails(run(z))
z <- x;z$samples <- z$samples[1,,drop=FALSE];z$samples[["Protein column"]] <- "Intensity";fails(run(z))
fails(run(x,protein_quantity="LFQ intensity"));fails(run(x,protein_quantity="automatic"))
fails(run(x,exclude_only_identified_by_site="true"))
z <- x;z$proteins$id[2] <- "0";fails(run(z))
z <- x;names(z$proteins)[names(z$proteins)=="Protein IDs"] <- "Other IDs";fails(run(z))

# Optional flags distinguish absent columns from present unflagged cells; there are no name heuristics.
z <- x;z$peptides$Reverse <- z$peptides[["Potential contaminant"]] <- NULL
z$proteins$Reverse <- z$proteins[["Potential contaminant"]] <- z$proteins[["Only identified by site"]] <- NULL
z$peptides[["Mod. peptide IDs"]] <- NULL
z$proteins[["Protein IDs"]][1] <- "CON__P1;REV__P2"
d <- run(z)
stopifnot(all(is.na(d$observations$contaminant_flag)),all(is.na(d$observations$decoy_flag)),
  all(is.na(d$entities$contaminant_flag)),all(is.na(d$entities$decoy_flag)),
  all(is.na(d$entities$only_identified_by_site_flag)),!any(d$observations$excluded),!any(d$entities$excluded),
  !any(d$observations$modifications_available),identical(d$source_peptides,z$peptides))
z <- x;z$peptides$Reverse[1] <- NA;stopifnot(!run(z)$observations$decoy_flag[1])
for (column in c("Potential contaminant","Reverse","Only identified by site")) {
  z <- x;z$proteins[[column]][1] <- "False";fails(run(z))
}
z <- x;z$peptides$Reverse[1] <- "-";fails(run(z))
z <- x;z$peptides$Sequence[7] <- "PEP[+16]TIDE";d <- run(z)
stopifnot(d$observations$exclusion_reason[7]=="unsupported_sequence_format;source_contaminant_flag")
cat("MaxQuant input tests passed.\n")
