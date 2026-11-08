# Generate one browsable, entirely synthetic report using the normal analysis stages.
# Rscript examples/linked_report/generate.R NEW_DIRECTORY
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L || file.exists(args[1])) stop("Supply a new output directory.")
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(script),"../.."),winslash="/",mustWork=TRUE)
for(module in c("mapping","cleavage","report","input_adapter","sample_ratios","normalization_io","statistics","statistics_io"))
  source(file.path(root,"R",paste0(module,".R")),encoding="UTF-8")
dir.create(args[1],recursive=TRUE)
out <- normalizePath(args[1],winslash="/")
inputs <- file.path(out,"inputs")
dir.create(inputs)
set.seed(751)
conditions <- rep(c("Control","Treatment"),each=6)
samples <- paste0(conditions,"_",rep(1:6,2))
proteins <- data.frame(Group=c("G1","G2"),Members=c("DEMO_A","DEMO_B"))
peptides <- data.frame(
  Observation=paste0("O",1:10),
  Sequence=c("MAAADFY","LQSMVR","MAAADFYLQSMVR","LQSMVR","GITPEPTK","LEFGK","MNPQTVK","LLDVAK","VTTASEK","AFY"),
  Modifications=c("","","","Oxidation [M4]",rep("",6)),
  Group=rep(c("G1","G2"),c(6,4))
)
effects <- c(1.0,1.4,-0.6,0.7,0,0,0,0,0,0.4)
for(i in seq_along(samples)){
  abundance <- c(25,24) + rnorm(2,0,0.2)
  proteins[[samples[i]]] <- abundance
  peptides[[samples[i]]] <- abundance[match(peptides$Group,proteins$Group)] - 5 +
    effects*(conditions[i]=="Treatment") + rnorm(nrow(peptides),0,0.24)
}
peptides$Control_4[2] <- NA
peptide_path <- file.path(inputs,"peptides.csv")
protein_path <- file.path(inputs,"proteins.csv")
pm_write_csv(peptides,peptide_path)
pm_write_csv(proteins,protein_path)
p <- pm_read_table(peptide_path)
q <- pm_read_table(protein_path)
fasta <- file.path(inputs,"synthetic.fasta")
writeLines(c(">DEMO_A Synthetic alpha GN=ALPHA","MAAADFYLQSMVRGITPEPTKAALEFGKLQSMVR",
             ">DEMO_B Synthetic beta GN=BETA","MNPQTVKLLDVAKVTTASEKAFY"),fasta)
mapped <- pm_map(p,pm_read_fasta(fasta),scope="all",sequence_column="Sequence")
mapped <- pm_build_cleavage(mapped,list(enzyme="trypsin",lysine_blocked=FALSE),digest::digest(file=fasta,algo="sha256"))
stopifnot(any(with(mapped$cleavage$sites,accession=="DEMO_A" & bond_left==7 & neo_n_count==2 & neo_c_count==1 & spanning_count==1)))
label <- "Synthetic demonstration (6 control and 6 treatment samples)"
pm_export(mapped,file.path(inputs,"mapping"),file.path(root,"inst","templates","report_template.html"),
          list(peptides=peptide_path,fasta=fasta),metadata=list(dataset_label=label),excel=FALSE)
metadata <- data.frame(Sample=samples,Condition=conditions,Use="Yes","Peptide column"=samples,"Protein column"=samples,check.names=FALSE)
spec <- list(observation_level="peptide_observation",link_separator=NULL,
             peptide=list(id="Observation",sequence="Sequence",modifications="Modifications",protein_links="Group",scale="log2",quantity_type="Synthetic log2 intensity"),
             protein=list(id="Group",accessions="Members",entity_type="protein_group",scale="log2",quantity_type="Synthetic log2 intensity"))
data <- pm_import_quant_tables(p,q,spec,metadata)
ratios <- pm_sample_ratios(data)
norm <- file.path(inputs,"ratios")
pm_export_sample_ratios(data,ratios,norm,list(peptides=peptide_path,proteins=protein_path),excel=FALSE,notes=label)
units <- data.frame(Sample=samples,Condition=conditions,Use="Yes","Experimental unit"=samples,check.names=FALSE)
contrasts <- data.frame(Name="Treatment_vs_Control",Treatment="Treatment",Control="Control")
results <- pm_test_ratios(ratios$ratios,units,contrasts)
pm_export_statistics(results,file.path(inputs,"statistics"),
  list(ratios=file.path(norm,"sample_ratios.csv"),ratio_manifest=file.path(norm,"normalization_manifest.json")),excel=FALSE,notes=label)
config <- file.path(inputs,"report.json")
jsonlite::write_json(list(mapping_dir="mapping",ratios_dir="ratios",statistics_dir="statistics",output_dir="../report"),config,auto_unbox=TRUE)
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
status <- system2(rscript,c(shQuote(file.path(root,"scripts","refresh_report.R")),"--config",shQuote(config)))
if(status!=0L) stop("Report generation failed.")
