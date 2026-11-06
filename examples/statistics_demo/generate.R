# Reproducible known-truth example; no experimental data or inference about biology.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: Rscript examples/statistics_demo/generate.R NEW_OUTPUT_DIRECTORY")
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),winslash="/",mustWork=TRUE)
root <- dirname(dirname(dirname(script)))
source(file.path(root,"R","report.R"),encoding="UTF-8")
out <- args[1]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE))) stop("Choose an empty example directory.")
dir.create(out,recursive=TRUE,showWarnings=FALSE);out <- normalizePath(out,winslash="/",mustWork=TRUE)
set.seed(7102026)
n <- 600L;ng <- 200L;ns <- 12L
samples <- c(paste0("Control_",1:6),paste0("Treatment_",1:6))
group <- rep(seq_len(ng),each=3L)
condition <- rep(c("Control","Treatment"),each=6L)
protein_shift <- c(rep(1,50),rep(0,150))
protein <- matrix(rep(rnorm(ng,25,1),ns),ng,ns)+matrix(rnorm(ng*ns,0,.3),ng,ns)+outer(protein_shift,as.numeric(condition=="Treatment"))
effect <- rep(0,n);effect[301:360] <- 1.5;effect[361:420] <- -1.5
sd <- runif(n,.2,.6)
noise <- matrix(rnorm(n*ns),n,ns)*sd
shared <- matrix(rnorm(ng*ns,0,.12),ng,ns)
ratio <- matrix(rep(rnorm(n,-5,1),ns),n,ns)+noise+shared[group,]+outer(effect,as.numeric(condition=="Treatment"))
peptide <- protein[group,]+ratio
# MCAR missingness is explicitly synthetic; no imputation is performed downstream.
peptide[matrix(runif(n*ns)<.04,n,ns)] <- NA_real_
protein[matrix(runif(ng*ns)<.01,ng,ns)] <- NA_real_
alphabet <- strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
sequence <- replicate(n,paste(sample(alphabet,14,replace=TRUE),collapse=""))
sequence[1:2] <- "PEPTIDEMK"
peptides <- data.frame(Observation=sprintf("P%04d",seq_len(n)),Sequence=sequence,
  Group=sprintf("G%03d",group),Modifications=c("","Oxidation [M8]",rep("",n-2)),Contaminant=FALSE,check.names=FALSE)
for(i in seq_len(ns)) peptides[[samples[i]]] <- peptide[,i]
proteins <- data.frame(Group=sprintf("G%03d",seq_len(ng)),Members=sprintf("SYNTHETIC_%03d",seq_len(ng)),check.names=FALSE)
for(i in seq_len(ns)) proteins[[samples[i]]] <- protein[,i]
constant <- proteins[1,,drop=FALSE];constant$Group <- "G201";constant$Members <- "SYNTHETIC_CONSTANT"
constant[,samples] <- 20;proteins <- rbind(proteins,constant)
diagnostics <- peptides[1:6,,drop=FALSE]
diagnostics$Observation <- c("ambiguous","low_support","constant","contaminant","no_link","unknown_link")
diagnostics$Group <- c("G001;G002","G003","G201","G004","","G999")
diagnostics$Sequence <- c("PEPTIDEK","PEPTIDER","ACDEFGHK","ACDEFGHR","KLMNPQR","KLMNPQK")
diagnostics$Modifications <- "";diagnostics$Contaminant <- c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE)
diagnostics[2,samples[1:5]] <- NA_real_
diagnostics[3,samples] <- c(rep(15,6),rep(17,6))
peptides <- rbind(peptides,diagnostics)
truth <- data.frame(observation_id=peptides$Observation,true_log2_effect=c(effect,rep(NA_real_,6)),
  scenario=c(ifelse(effect!=0,"peptide_change",ifelse(protein_shift[group]!=0,"protein_only_change","null")),diagnostics$Observation))
metadata <- data.frame(Sample=samples,Condition=condition,Use="Yes","Peptide column"=samples,"Protein column"=samples,check.names=FALSE)
stats_samples <- data.frame(Sample=samples,Condition=condition,Use="Yes","Experimental unit"=paste0("Unit_",seq_len(ns)),check.names=FALSE)
comparisons <- data.frame(Name="Treatment_vs_Control",Treatment="Treatment",Control="Control")
for(name in c("peptides","proteins","truth","metadata","stats_samples","comparisons")) pm_write_csv(get(name),file.path(out,paste0(name,".csv")))
normalization <- list(peptides="peptides.csv",proteins="proteins.csv",samples="metadata.csv",output_dir="ratios",excel=FALSE,
  spec=list(observation_level="peptide_observation",link_separator=";",
    peptide=list(id="Observation",sequence="Sequence",protein_links="Group",modifications="Modifications",contaminant="Contaminant",scale="log2",quantity_type="Synthetic log2 peptide quantity"),
    protein=list(id="Group",accessions="Members",entity_type="protein_group",scale="log2",quantity_type="Synthetic log2 group quantity")),
  notes="Synthetic seed 7102026; known effects and MCAR missingness, not experimental biology.")
statistics <- list(ratios_dir="ratios",samples="stats_samples.csv",comparisons="comparisons.csv",output_dir="statistics",
  design="independent_units",excel=TRUE,notes=c("Synthetic seed 7102026: 600 modeled observations plus six diagnostic examples, 6 independent units per condition.",
  "Protein-only changes have true adjusted effect zero; 60 positive and 60 negative peptide-specific effects have magnitude 1.5 log2 units. Truth is separate from fitted input.",
  "Shared protein variation and group-specific ratio noise induce dependence. MCAR missingness is generated explicitly; this example is not a biological benchmark."))
jsonlite::write_json(normalization,file.path(out,"normalize.json"),auto_unbox=TRUE,pretty=TRUE)
jsonlite::write_json(statistics,file.path(out,"statistics.json"),auto_unbox=TRUE,pretty=TRUE)
cat("Synthetic inputs and explicit configurations written to",out,"\n")
