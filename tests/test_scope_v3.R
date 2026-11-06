script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-dirname(dirname(normalizePath(script,winslash="/")))
for(name in c("mapping","cleavage","pair_statistics","pair_workflow"))source(file.path(root,"R",paste0(name,".R")),encoding="UTF-8")
td<-tempfile("scope_v3_");dir.create(td)
fa<-file.path(td,"reference.fasta")
writeLines(c(">sp|P1|ONE","MAARACDEFGHIKLMNRSTVWYK",">sp|P2|ALTERNATIVE","MAARACDEFGHIKLMNRSTVWYK"),fa)
p<-data.frame(Sequence=c("CDEFGHIK","DEFGHIK","ACDEFGHIK","CDEFGHI","AAR","MAAR"),
  Modifications="","Protein Accessions"="P1;P2",Chosen="P1",check.names=FALSE)
quant<-rbind(c(8,10,11,20,21,24),c(4,6,5,15,12,14),c(100,110,105,108,99,112),
  c(8,9,8,9,8,9),c(10,11,9,15,13,16),c(20,22,19,21,24,23))
for(i in 1:6)p[[paste0("Abundance: S",i)]]<-quant[,i]
x<-pm_map(p,pm_read_fasta(fa),sequence_column="Sequence")
legacy<-pm_build_cleavage(x,list(enzyme="trypsin",lysine_blocked=FALSE),digest::digest(file=fa,algo="sha256"))
v3<-pm_build_cleavage(x,list(enzyme="trypsin",lysine_blocked=FALSE,analysis_scope="specificity_v3"),digest::digest(file=fa,algo="sha256"))
stopifnot(identical(x$measurements,v3$measurements),identical(x$matches,v3$matches),
  is.null(legacy$cleavage$placements),nrow(v3$cleavage$termini)==nrow(legacy$cleavage$termini))
audit<-v3$cleavage$placements
stopifnot(all(audit$specificity[audit$peptide_id=="pep_0000001"]=="semi-specific"),
  all(audit$specific_terminus[audit$peptide_id=="pep_0000001"]=="C"),
  all(audit$specificity[audit$peptide_id=="pep_0000003"]=="specific"),
  all(audit$specificity[audit$peptide_id=="pep_0000004"]=="unspecific"),
  all(audit$specificity[audit$peptide_id=="pep_0000005"]=="specific"),
  all(audit$n_processing_context[audit$peptide_id=="pep_0000005"]=="possible_initiator_met_removal"),
  !"pep_0000004" %in% v3$cleavage$evidence$peptide_id,
  all(v3$cleavage$excluded_placements$peptide_id=="pep_0000004"))
# N labels do not redefine analytical sequence specificity. Named unknown profiles
# remain unknown and cannot be silently converted to a tryptic rule.
labelp<-p;labelp$Modifications[4]<-"TMTpro [N-Term]"
labelled<-pm_build_cleavage(pm_map(labelp,pm_read_fasta(fa),sequence_column="Sequence"),
  list(enzyme="trypsin",lysine_blocked=FALSE,pre_digest_n_label="TMTpro",analysis_scope="specificity_v3"),digest::digest(file=fa,algo="sha256"))
stopifnot(!"pep_0000004" %in% labelled$cleavage$evidence$peptide_id,
  any(labelled$cleavage$termini$is_candidate_before_scope[labelled$cleavage$termini$peptide_id=="pep_0000004"]))
unknown<-pm_build_cleavage(x,list(enzyme="unknown",analysis_scope="specificity_v3"),digest::digest(file=fa,algo="sha256"))
stopifnot(all(unknown$cleavage$placements$specificity=="unknown"),nrow(unknown$cleavage$excluded_placements)==0)
samples<-data.frame(Sample=paste0("S",1:6),Condition=rep(c("Control","Treatment"),each=3),Use="Yes",
  "Peptide column"=paste0("Abundance: S",1:6),"Experimental unit"=paste0("Animal",1:6),check.names=FALSE)
design<-list(samples=samples,comparisons=data.frame(Name="T-C",Treatment="Treatment",Control="Control"))
unique_scope<-pm_pair_workflow(v3,design,assignment_mode="unique",test=FALSE)
stopifnot(!any(unique_scope$links$scope_eligible),any(grepl("Ambiguous protein assignment",unique_scope$excluded$reason,fixed=TRUE)))
representative<-pm_pair_workflow(v3,design,representative_column="Chosen",assignment_mode="representative")
stopifnot(any(representative$measurement_results$tested),
  all(representative$results$accession=="P1"),
  !any(representative$results$bond_left==1L),
  any(grepl("initiator-Met",representative$excluded$reason,fixed=TRUE)),
  any(grepl("Outside shared master",representative$excluded$reason,fixed=TRUE)),
  !anyDuplicated(representative$measurement_results[,c("measurement_pair_id","family","comparison")]))
# Repeating annotation links must not change measurements, moderation or BH.
dup<-v3;dup$cleavage$evidence<-rbind(dup$cleavage$evidence,dup$cleavage$evidence)
again<-pm_pair_workflow(dup,design,"Chosen","representative")
stopifnot(identical(representative$measurement_results,again$measurement_results))
# Position 1|2 alone is not initiator-Met processing: require an actual initial M.
fa_nonmet<-file.path(td,"nonmet.fasta")
writeLines(c(">sp|P1|ONE","AAARACDEFGHIKLMNRSTVWYK",">sp|P2|ALTERNATIVE","AAARACDEFGHIKLMNRSTVWYK"),fa_nonmet)
nonmet_p<-p;nonmet_p$Sequence[6]<-"AAAR"
nonmet<-pm_build_cleavage(pm_map(nonmet_p,pm_read_fasta(fa_nonmet),sequence_column="Sequence"),
  list(enzyme="trypsin",lysine_blocked=FALSE,analysis_scope="specificity_v3"),digest::digest(file=fa_nonmet,algo="sha256"))
nonmet_out<-pm_pair_workflow(nonmet,design,"Chosen","representative",test=FALSE)
stopifnot(any(nonmet_out$links$scope_eligible & nonmet_out$links$bond_left==1L))
# Missing modification metadata is distinct from an explicitly unlabelled form.
state_p<-p;state_p$Modifications[1]<-NA_character_
state_map<-pm_build_cleavage(pm_map(state_p,pm_read_fasta(fa),sequence_column="Sequence"),
  list(enzyme="trypsin",lysine_blocked=FALSE,pre_digest_n_label="TMTpro",analysis_scope="specificity_v3"),digest::digest(file=fa,algo="sha256"))
state_out<-pm_pair_workflow(state_map,design,"Chosen","representative",test=FALSE)
stopifnot(all(state_out$links$family[state_out$links$neo_peptide_id=="pep_0000001"]=="N_unresolved_exploratory"),
  all(state_out$links$family[state_out$links$neo_peptide_id=="pep_0000002"]=="N_unlabelled_exploratory"))
# Master assignments are sets. Their intersection defines eligible attributions,
# while one directed observation pair remains one test per family/contrast.
bad<-v3;bad$source_peptides$Chosen[1]<-"P2"
conflict<-pm_pair_workflow(bad,design,"Chosen","representative",test=FALSE)
stopifnot(any(grepl("No shared master",conflict$excluded$reason)))
bad$source_peptides$Chosen[1]<-"P1;P2"
multiple<-pm_pair_workflow(bad,design,"Chosen","representative",test=FALSE)
stopifnot(any(multiple$links$scope_eligible & multiple$links$neo_peptide_id=="pep_0000001"),
 all(multiple$links$common_master_accessions[multiple$links$neo_peptide_id=="pep_0000001"]=="P1"))
# Two list-valued assignments may share only P1, irrespective of input order.
intersection_one<-v3;intersection_one$source_peptides$Chosen<-"P1;P3"
intersection_one$source_peptides$Chosen[1:2]<-"P2;P1"
one<-pm_pair_workflow(intersection_one,design,"Chosen","representative")
stopifnot(all(one$results$accession[one$results$spanning_peptide_id=="pep_0000003"]=="P1"),
 identical(one$measurement_results,representative$measurement_results))
intersection_two<-v3;intersection_two$source_peptides$Chosen<-"P2;P1;P1"
two<-pm_pair_workflow(intersection_two,design,"Chosen","representative")
stopifnot(setequal(two$results$accession,c("P1","P2")),
 all(two$results$common_master_accessions=="P1;P2"),all(two$results$common_master_count==2L),
 all(two$results$attribution_ambiguous),all(two$results$eligible_mapping_count==2L),
 identical(two$measurement_results,representative$measurement_results),
 !anyDuplicated(two$measurement_results[,c("measurement_pair_id","family","comparison")]))
two_dup<-intersection_two;two_dup$cleavage$evidence<-rbind(two_dup$cleavage$evidence,two_dup$cleavage$evidence)
stopifnot(identical(pm_pair_workflow(two_dup,design,"Chosen","representative")$measurement_results,two$measurement_results))
# First sequence occurrence is deterministic within one protein record; later
# occurrences remain in the audit and cannot rescue an unsupported first pair.
repeated_fa<-file.path(td,"repeated.fasta")
writeLines(c(">sp|P1|ONE",paste0(rep("MAARACDEFGHIKLMNRSTVWYK",2),collapse="")),repeated_fa)
repeat_p<-p[1:3,];repeat_p[["Protein Accessions"]]<-"P1"
repeat_map<-pm_build_cleavage(pm_map(repeat_p,pm_read_fasta(repeated_fa),sequence_column="Sequence"),
 list(enzyme="trypsin",lysine_blocked=FALSE,analysis_scope="specificity_v3"),digest::digest(file=repeated_fa,algo="sha256"))
repeat_out<-pm_pair_workflow(repeat_map,design,"Chosen","representative")
stopifnot(any(repeat_out$links$scope_eligible),any(!repeat_out$links$assignment_eligible),
 all(repeat_out$links$attribution_ambiguous),nrow(repeat_out$measurement_results)>0L,
 all(repeat_out$links$neo_placement_count==2L),all(repeat_out$links$spanning_placement_count==2L),
 all(repeat_out$links$bond_left[repeat_out$links$scope_eligible]<22L),
 all(grepl("First sequence occurrences",repeat_out$excluded$reason)),
 !anyDuplicated(repeat_out$measurement_results[,c("measurement_pair_id","family","comparison")]))
repeat_unique<-pm_pair_workflow(repeat_map,design,assignment_mode="unique",test=FALSE)
stopifnot(identical(repeat_unique$links$assignment_eligible,repeat_out$links$assignment_eligible))
shuffled<-repeat_map;shuffled$matches<-shuffled$matches[rev(seq_len(nrow(shuffled$matches))),,drop=FALSE]
shuffled_out<-pm_pair_workflow(shuffled,design,"Chosen","representative")
stopifnot(identical(repeat_out$links,shuffled_out$links),
 identical(repeat_out$measurement_results,shuffled_out$measurement_results))
# Two FASTA records under one accession are not sequence-occurrence ties.
duplicate_record<-repeat_map
extra<-duplicate_record$matches;extra$record_id<-max(extra$record_id)+1L
duplicate_record$matches<-rbind(duplicate_record$matches,extra)
duplicate_out<-pm_pair_workflow(duplicate_record,design,"Chosen","representative",test=FALSE)
stopifnot(!any(duplicate_out$links$assignment_eligible),
 all(grepl("Ambiguous placement",duplicate_out$excluded$reason)))
disjoint<-v3;disjoint$source_peptides$Chosen[1:2]<-"P2;P3"
disjoint$source_peptides$Chosen[3]<-"P1;P4"
disjoint_out<-pm_pair_workflow(disjoint,design,"Chosen","representative",test=FALSE)
stopifnot(!any(disjoint_out$links$scope_eligible[disjoint_out$links$spanning_peptide_id=="pep_0000003"]))
fails<-function(expr)stopifnot(inherits(try(expr,silent=TRUE),"try-error"))
state_override<-data.frame(peptide_id="pep_0000002",state="labelled")
fails(pm_pair_workflow(state_map,design,"Chosen","representative",test=FALSE,n_terminal_states=state_override))
acetyl_p<-p;acetyl_p$Modifications[1]<-"1xAcetyl [N-Term]"
acetyl_map<-pm_build_cleavage(pm_map(acetyl_p,pm_read_fasta(fa),sequence_column="Sequence"),
  list(enzyme="trypsin",lysine_blocked=FALSE,pre_digest_n_label="TMTpro",analysis_scope="specificity_v3"),digest::digest(file=fa,algo="sha256"))
state_override$peptide_id<-"pep_0000001"
fails(pm_pair_workflow(acetyl_map,design,"Chosen","representative",test=FALSE,n_terminal_states=state_override))
mass_p<-acetyl_p;mass_p$Modifications[1]<-"Mass shift 304.2071 [N-Term]"
mass_map<-pm_build_cleavage(pm_map(mass_p,pm_read_fasta(fa),sequence_column="Sequence"),
  list(enzyme="trypsin",lysine_blocked=FALSE,pre_digest_n_label="TMTpro",analysis_scope="specificity_v3"),digest::digest(file=fa,algo="sha256"))
resolved<-pm_pair_workflow(mass_map,design,"Chosen","representative",test=FALSE,n_terminal_states=state_override)
stopifnot(all(resolved$links$family[resolved$links$neo_peptide_id=="pep_0000001"]=="N_TMT_primary"))
fails(pm_pair_workflow(v3,design,"not present","representative"))
fails(pm_pair_workflow(v3,design,"Chosen","representative",abundance_type="normalized"))
bad_design<-design;bad_design$samples$Use[2]<-"unknown"
fails(pm_pair_workflow(v3,bad_design,"Chosen","representative"))
bad_design<-design;bad_design$samples[["Experimental unit"]][2]<-"Animal1"
fails(pm_pair_workflow(v3,bad_design,"Chosen","representative"))
bad_design<-design;bad_design$comparisons<-rbind(design$comparisons,design$comparisons)
bad_design$comparisons$Name[2]<-"Same contrast with another name"
fails(pm_pair_workflow(v3,bad_design,"Chosen","representative"))
consistent_units<-samples[,c("Sample","Experimental unit","Condition","Use")]
consistent<-pm_pair_workflow(v3,design,"Chosen","representative",experimental_units=consistent_units[6:1,])
stopifnot(identical(representative$measurement_results,consistent$measurement_results))
bad_units<-consistent_units;bad_units$Condition[2]<-"WRONG"
fails(pm_pair_workflow(v3,design,"Chosen","representative",experimental_units=bad_units))
bad_units<-consistent_units;bad_units$Use[2]<-"No"
fails(pm_pair_workflow(v3,design,"Chosen","representative",experimental_units=bad_units))
bad_units$Use[2]<-NA_character_
fails(pm_pair_workflow(v3,design,"Chosen","representative",experimental_units=bad_units))
# Missing/zero measurements remain unavailable and cannot manufacture a finite FC.
missing<-v3;selected<-missing$measurements$peptide_id=="pep_0000001" & missing$measurements$source_column %in% samples[["Peptide column"]][1:3]
missing$measurements$value[selected]<-NA_real_;missing$measurements$value_status[selected]<-"missing"
missing_result<-pm_pair_workflow(missing,design,"Chosen","representative")
z<-missing_result$contrasts[missing_result$contrasts$neo_peptide_id=="pep_0000001",]
stopifnot(all(is.na(z$log2_fc)),all(z$status=="neo_quantified_only_in_treatment"))
cat("PASS v3 specificity, unspecific exclusion, natural-boundary audit, representative scope, measurement deduplication and missingness.\n")
