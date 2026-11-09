# Descriptive spanning ends, shared by workbook and HTML reporting.
pm_report_spanning_ends <- function(payload) {
  out <- data.frame(site_id=character(),spanning_peptide_id=character(),
    ends=character(),terminal=character(),position=character(),stringsAsFactors=FALSE)
  evidence <- payload$cleavage$evidence; terms <- payload$cleavage$termini
  if(!is.data.frame(evidence)||!nrow(evidence)) return(out)
  evidence <- unique(evidence[evidence$role %in% "spanning",
    intersect(c("site_id","peptide_id","record_id","accession","start","end"),names(evidence)),drop=FALSE])
  if(!nrow(evidence))return(out)
  key <- function(x,cols) do.call(paste,c(x[cols],sep="\034"))
  index <- function(x,cols) split(seq_len(nrow(x)),key(x,cols))
  end_labels <- character();term_index<-list()
  required<-c("record_id","peptide_id","start","end","role")
  if(is.data.frame(terms)&&nrow(terms)&&all(required %in% names(terms))){
    value<-function(name,default) if(name %in% names(terms))terms[[name]] else rep(default,nrow(terms))
    ps<-as.character(value("processing_status",""));cl<-as.character(value("classification",""))
    labelled<-value("pre_digest_n_label",NA);compatible<-value("digestion_compatible",NA)
    end_labels<-ifelse(ps %in% "reference_end"|cl %in% "reference_end","protein end",
      ifelse(ps %in% "possible_initiator_met_removal","possible Met removal",
        ifelse(terms$role %in% "N" & labelled %in% TRUE,"pre-digestion N label",
          ifelse(compatible %in% TRUE,"digestion-compatible",
            ifelse(compatible %in% FALSE,"not digestion-compatible","unknown")))))
    term_index<-index(terms,required)
  }
  groups<-index(evidence,c("site_id","peptide_id"))
  rows<-lapply(groups,function(ii){
    descriptions<-character();labels<-character()
    for(j in ii){
      side<-vapply(c("N","C"),function(role){
        kk<-paste(evidence$record_id[j],evidence$peptide_id[j],evidence$start[j],evidence$end[j],role,sep="\034")
        tt<-term_index[[kk]]
        if(length(tt))paste(unique(end_labels[tt]),collapse="/") else "unknown"
      },character(1))
      labels<-c(labels,unlist(strsplit(side,"/",fixed=TRUE)))
      text<-paste0("N: ",side[1],"; C: ",side[2])
      if(length(ii)>1L){
        acc<-if("accession" %in% names(evidence))evidence$accession[j] else evidence$record_id[j]
        text<-paste0(acc," [",evidence$record_id[j]," ",evidence$start[j],"-",evidence$end[j],"]: ",text)
      }
      descriptions<-c(descriptions,text)
    }
    flag<-if(any(labels %in% c("pre-digestion N label","not digestion-compatible")))"Yes" else if(!length(labels)||any(labels %in% "unknown"))"Unknown" else "No"
    data.frame(site_id=evidence$site_id[ii[1]],spanning_peptide_id=evidence$peptide_id[ii[1]],
      ends=paste(unique(descriptions),collapse=" | "),terminal=flag,
      position=paste(unique(paste(evidence$start[ii],evidence$end[ii],sep="-")),collapse="; "),stringsAsFactors=FALSE)
  })
  out<-do.call(rbind,rows);rownames(out)<-NULL
  out
}

# Reader-facing tables, derived from unchanged observations and pair tests.
# This mirrors the frozen report representative rule; no models are fitted here.
pm_export_tables <- function(payload) {
  empty <- function() data.frame(Note=character(),check.names=FALSE)
  bind <- function(xs) {
    xs<-Filter(function(x)is.data.frame(x)&&nrow(x)>0,xs)
    if(!length(xs)) return(empty())
    cols<-unique(unlist(lapply(xs,names)));xs<-lapply(xs,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-NA;x[cols]})
    do.call(rbind,xs)
  }
  key <- function(x,cols) {if(!is.data.frame(x)||!nrow(x))return(character());do.call(paste,c(x[cols],sep="\034"))}
  index <- function(x,cols) split(seq_len(nrow(x)),key(x,cols))
  collapse <- function(x) paste(unique(as.character(x[!is.na(x)&nzchar(as.character(x))])),collapse="; ")
  mapping_only<-is.null(payload$cleavage)
  pv<-if(mapping_only&&!is.null(payload$sample_quantities))payload$sample_quantities else payload$pair_view
  ratios<-pv$ratios;sites<-payload$cleavage$sites;evidence<-payload$cleavage$evidence;terms<-payload$cleavage$termini
  source_column<-function(column){
    q<-payload$column_mapping_audit$quantities
    if(is.data.frame(q)&&all(c("canonical_column","source_column") %in% names(q))){
      j<-match(column,q$canonical_column);if(!is.na(j))return(q$source_column[j])
    }
    column
  }
  if(is.null(evidence))evidence<-data.frame()
  if(is.null(terms))terms<-data.frame()
  paircols<-c("site_id","role","neo_peptide_id","spanning_peptide_id")
  prepared<-pm_report_pair_summaries(payload)
  eligible_rows<-if(nrow(prepared))which(prepared$assignment_eligible) else integer()
  all<-if(length(eligible_rows))prepared[eligible_rows,,drop=FALSE] else empty()
  chosen<-if(length(eligible_rows))all[all$Representative,,drop=FALSE] else empty()
  # One canonical grouped table also backs filtered downloads. It contains every
  # representative; the public default workbook applies only the display filter.
  public_results<-function(x){
    if(!"Comparison" %in% names(x))return(x)
    comps<-unique(x$Comparison)
    if(length(comps)==1L && !is.na(comps) && nzchar(comps)){
      for(side in c("Control","Treatment")){
        label<-unique(x[[side]]);total<-unique(x[[paste(side,"total")]])
        if(length(label)==1L&&length(total)==1L)
          names(x)[names(x)==paste(side,"measured")]<-paste0(label," measured (of ",total,")")
      }
    }
    x[setdiff(names(x),c("Control","Control total","Treatment","Treatment total","Site IDs","Export key","Passes default coverage"))]
  }
  # Ends of the spanning form at its placement(s) for the focal bond. This is
  # descriptive: representative selection, ratios and tests are unchanged. A
  # spanning form is reported as a terminal form when one of its own ends was not
  # produced by the declared digestion (a non-compatible end, or the declared
  # pre-digestion N-terminal label away from a protein end).
  reporting<-pv$spanning_end_reporting
  if(!is.data.frame(reporting))reporting<-pm_report_spanning_ends(payload)
  report_keys<-key(reporting,c("site_id","spanning_peptide_id"))
  terminal_form<-function(labels){
    if(!length(labels)||any(labels=="unknown"))return(if(any(labels %in% c("pre-digestion N label","not digestion-compatible")))"Yes" else "Unknown")
    if(any(labels %in% c("pre-digestion N label","not digestion-compatible")))"Yes" else "No"
  }
  spanning_ends<-function(site,peptide){
    at<-match(paste(site,peptide,sep="\034"),report_keys)
    ends<-reporting$ends[at];terminal<-reporting$terminal[at]
    ends[is.na(ends)]<-"Not available";terminal[is.na(terminal)]<-"Unknown"
    list(ends=ends,terminal=terminal)
  }
  all_index<-if("site_id" %in% names(all))index(all,c("site_id","role","comparison")) else list()
  position_index<-if(nrow(evidence))index(evidence,c("site_id","role","peptide_id")) else list()
  position_values<-vapply(position_index,function(ii)collapse(paste(evidence$start[ii],evidence$end[ii],sep="-")),character(1))
  display<-function(x,group=FALSE){
    if(!"site_id" %in% names(x))return(empty())
    stats<-pv$statistics;sk<-key(stats,c(paircols,"comparison"));sx<-split(seq_along(sk),sk);refs<-payload$site_summary
    x$`P-value`<-NA_real_;x$`BH-adjusted P`<-NA_real_;x$`Test family`<-"";x$`Test status`<-"Not tested"
    reasons<-c(insufficient_replicates="Insufficient replicates",zero_residual_variance="No residual variation",insufficient_features_for_moderation="Insufficient features for testing",numeric_failure="Numeric failure",model_numeric_failure="Model numeric failure")
    excluded<-pv$excluded_statistics;excluded_keys<-key(excluded,c(paircols,"comparison"))
    xkeys<-key(x,c(paircols,"comparison"));matched<-match(xkeys,names(sx));counts<-lengths(sx)[matched];counts[is.na(counts)]<-0L
    single<-which(counts==1L)
    if(length(single)){
      ii<-vapply(sx[matched[single]],function(z)z[1],integer(1))
      x$`Test family`[single]<-stats$family[ii];status<-stats$status[ii]
      tested<-status %in% "tested";rows<-single[tested]
      x$`P-value`[rows]<-stats$p_value[ii[tested]];x$`BH-adjusted P`[rows]<-stats$adjusted_p_value[ii[tested]]
      x$`Test status`[rows]<-"Tested"
      statuses<-ifelse(status %in% names(reasons),reasons[status],gsub("_"," ",status,fixed=TRUE))
      x$`Test status`[single[!tested]]<-statuses[!tested]
    }
    x$`Test status`[counts>1L]<-"Multiple test families"
    absent<-which(counts==0L);j<-match(xkeys[absent],excluded_keys);present<-!is.na(j)
    if(any(present))x$`Test status`[absent[present]]<-excluded$reason[j[present]]
    x$Gene<-if(is.data.frame(payload$proteins)&&"gene_symbol" %in% names(payload$proteins))payload$proteins$gene_symbol[match(x$record_id,payload$proteins$record_id)] else ""
    x$`Neo sequence`<-payload$peptides$sequence[match(x$neo_peptide_id,payload$peptides$peptide_id)]
    x$`Spanning sequence`<-payload$peptides$sequence[match(x$spanning_peptide_id,payload$peptides$peptide_id)]
    x$MEROPS<-if(!is.null(refs))refs$merops_summary[match(x$site_id,refs$site_id)] else "Not supplied"
    x$UniProt<-if(!is.null(refs))refs$uniprot_summary[match(x$site_id,refs$site_id)] else "Not supplied"
    x$Bond<-paste(x$bond_left,x$bond_right,sep=" | ")
    x$`Cleavage window`<-if("context_display" %in% names(sites))sites$context_display[match(x$site_id,sites$site_id)] else NA_character_
    x$`Neo modifications`<-payload$peptides$modifications[match(x$neo_peptide_id,payload$peptides$peptide_id)]
    x$`Spanning modifications`<-payload$peptides$modifications[match(x$spanning_peptide_id,payload$peptides$peptide_id)]
    x$Pattern<-ifelse(is.finite(x$log2_fc),ifelse(x$log2_fc>=1,"Ratio increased",ifelse(x$log2_fc<= -1,"Ratio decreased","Below ratio threshold")),x$status)
    x$Pattern[grepl(" quantified only in ",x$status,fixed=TRUE)]<-"Quantified in one condition"
    x$Pattern[x$status=="Spanning quantified; neo not quantified"]<-"Spanning quantified \u00b7 neo not quantified"
    x$Pattern[x$status=="Neo quantified; spanning not quantified"]<-"Neo quantified \u00b7 spanning not quantified"
    group_keys<-key(x,c("site_id","role","comparison"))
    x$`Available pairs`<-unname(lengths(all_index)[match(group_keys,names(all_index))]);x$`Available pairs`[is.na(x$`Available pairs`)]<-0L
    positions<-function(role,id){
      kk<-paste(x$site_id,role,id,sep="\034")
      v<-position_values[match(kk,names(position_values))];v[is.na(v)]<-"";unname(v)
    }
    x$`Neo position`<-positions(x$role,x$neo_peptide_id)
    x$`Spanning position`<-positions("spanning",x$spanning_peptide_id)
    span_ends<-spanning_ends(x$site_id,x$spanning_peptide_id)
    x$`Spanning ends`<-span_ends$ends;x$`Spanning terminal form`<-span_ends$terminal
    if(group){
      x<-bind(lapply(split(seq_len(nrow(x)),key(x,c("comparison","role","neo_peptide_id","spanning_peptide_id"))),function(ii){z<-x[ii[1],,drop=FALSE];z$accession<-collapse(x$accession[ii]);z$Gene<-collapse(x$Gene[ii]);z$Bond<-collapse(paste(x$accession[ii],x$Bond[ii]));z$site_id<-collapse(x$site_id[ii]);z$MEROPS<-collapse(x$MEROPS[ii]);z$UniProt<-collapse(x$UniProt[ii]);
        descriptions<-x$`Spanning ends`[ii]
        if(length(ii)>1L)descriptions<-paste0(x$accession[ii]," [",x$record_id[ii]," ",x$`Spanning position`[ii],"]: ",descriptions)
        z$`Spanning ends`<-collapse(descriptions);z$`Spanning terminal form`<-terminal_form(unlist(lapply(x$`Spanning terminal form`[ii],function(v)switch(v,Yes="pre-digestion N label",No="digestion-compatible","unknown"))));
        windows<-unique(x$`Cleavage window`[ii]);z$`Cleavage window`<-if(length(windows)==1L)windows else collapse(paste(x$accession[ii],x$`Cleavage window`[ii]));
        alternatives<-unique(unlist(lapply(ii,function(i)all_index[[key(x[i,,drop=FALSE],c("site_id","role","comparison"))]])));
        z$`Available pairs`<-length(unique(key(all[alternatives,,drop=FALSE],c("neo_peptide_id","spanning_peptide_id"))));
        signatures<-key(x[ii,,drop=FALSE],c("Test family","P-value","BH-adjusted P","Test status"));if(length(unique(signatures))>1){z$`P-value`<-NA_real_;z$`BH-adjusted P`<-NA_real_;z$`Test family`<-"";z$`Test status`<-"Varies by mapping; see details"};z}))
    }
    x$`Export key`<-vapply(seq_len(nrow(x)),function(i)as.character(jsonlite::toJSON(unname(as.character(x[i,c("comparison","role","neo_peptide_id","spanning_peptide_id")])) ,auto_unbox=FALSE,na="null")),character(1))
    x$role<-paste0("Neo-",x$role)
    fields<-c(comparison="Comparison",accession="Protein",Gene="Gene",Bond="Bond",role="Role",`Cleavage window`="Cleavage window",
      `Neo sequence`="Neo sequence",`Neo modifications`="Neo modifications",`Spanning sequence`="Spanning sequence",`Spanning modifications`="Spanning modifications",
      `Spanning ends`="Spanning ends",`Spanning terminal form`="Spanning terminal form",
      neo_log2_fc="Neo log2FC",spanning_log2_fc="Spanning log2FC",log2_fc="Ratio log2FC",Pattern="Pattern",
      control_common="Control measured",treatment_common="Treatment measured",`Available pairs`="Available pairs",status="Status",
      `P-value`="P-value",`BH-adjusted P`="BH-adjusted P",`Test family`="Test family",`Test status`="Test status",MEROPS="MEROPS",UniProt="UniProt")
    if(!group)fields<-c(fields,`Neo position`="Neo position",`Spanning position`="Spanning position",Representative="Representative")
    x$`Annotation IDs`<-x$site_id
    fields<-c(fields,`Annotation IDs`="Annotation IDs",neo_peptide_id="Neo form ID",spanning_peptide_id="Spanning form ID",control="Control",control_total="Control total",treatment="Treatment",treatment_total="Treatment total",site_id="Site IDs")
    if(group)fields<-c(fields,`Export key`="Export key",coverage_pass="Passes default coverage")
    pm_excel_select(x,fields)
  }
  # A source form is never expanded by the number of mapped bonds.
  ev<-evidence[evidence$is_candidate %in% TRUE|evidence$role %in% "spanning",,drop=FALSE]
  ids<-unique(ev$peptide_id);peps<-payload$peptides[match(ids,payload$peptides$peptide_id),,drop=FALSE]
  if(is.null(payload$cleavage)){ids<-payload$peptides$peptide_id;peps<-payload$peptides}
  ex<-split(seq_len(nrow(ev)),ev$peptide_id);span_sites<-unique(ev$site_id[ev$role=="spanning"])
  used<-if("neo_peptide_id" %in% names(all))unique(c(all$neo_peptide_id,all$spanning_peptide_id)) else character()
  if(length(ids)&&nrow(ev)){
    peps$roles<-vapply(ids,function(id)collapse(ifelse(ev$role[ex[[id]]]=="spanning","Spanning",paste0("Neo-",ev$role[ex[[id]]]))),character(1))
    peps$locations<-vapply(ids,function(id){z<-ev[ex[[id]],,drop=FALSE];collapse(paste(z$accession,paste(z$start,z$end,sep="-")))},character(1))
    peps$usage<-vapply(ids,function(id){if(id %in% used)return("Used in eligible pair")
      z<-ev[ex[[id]],,drop=FALSE];z<-z[z$role!="spanning",,drop=FALSE]
      if(nrow(z)&&any(z$site_id %in% span_sites))"Spanning available; no eligible quantified pair" else if(nrow(z))"No spanning peptide" else "Supporting spanning form"},character(1))
  }
  peptides<-pm_excel_select(peps,c(peptide_id="Form ID",sequence="Peptide sequence",modifications="Modifications",roles="Roles",locations="Protein positions",usage="Usage",source_row="Source row"))
  peptide_columns<-list()
  add_peptide_column<-function(label,context,value){
    peptide_columns[[length(peptide_columns)+1L]]<<-list(label=label,context=context,value=value)
  }
  if(is.data.frame(pv$measurements)&&nrow(pv$measurements)&&is.data.frame(pv$comparisons)&&nrow(pv$comparisons)){
    mm<-pv$measurements;mi<-split(seq_len(nrow(mm)),mm$peptide_id)
    for(ci in seq_len(nrow(pv$comparisons))){cp<-pv$comparisons[ci,,drop=FALSE];ss<-pv$samples[pv$samples$Use=="Yes",,drop=FALSE]
      cols<-lapply(c(cp$Control,cp$Treatment),function(cond)ss[["Peptide column"]][ss$Condition==cond]);totals<-lengths(cols)
      values<-lapply(ids,function(id){m<-mm[mi[[id]],,drop=FALSE];lapply(cols,function(cc){v<-m$value[m$source_column %in% cc&m$value_status %in% "observed"&is.finite(m$value)&m$value>0];list(n=length(v),mean=if(length(v))mean(log2(v)) else NA_real_)})})
      add_peptide_column(paste(cp$Name,"peptide log2FC"),paste("comparison",cp$Name),vapply(values,function(v)v[[2]]$mean-v[[1]]$mean,numeric(1)))
      for(j in 1:2){lab<-c(cp$Control,cp$Treatment)[j];add_peptide_column(paste0(cp$Name," ",lab," measured (of ",totals[j],")"),paste0("comparison ",cp$Name,"; ",c("control","treatment")[j]," ",lab),vapply(values,function(v)v[[j]]$n,integer(1)))}
    }
  }
  # Export measured quantities exactly as supplied to this report. Values are
  # not rounded, log-transformed, imputed or screened by ratio eligibility here.
  mm<-pv$measurements
  if(is.data.frame(mm)&&nrow(mm)){
    if(anyDuplicated(key(mm,c("peptide_id","source_column"))))stop("Duplicate form/sample quantities in workbook input.")
    sample_columns<-if(is.data.frame(pv$samples)&&"Peptide column" %in% names(pv$samples))unique(pv$samples[["Peptide column"]]) else unique(mm$source_column)
    sample_labels<-vapply(sample_columns,function(column){
      if(is.data.frame(pv$samples)&&"Sample" %in% names(pv$samples)){
        j<-match(column,pv$samples[["Peptide column"]]);s<-pv$samples$Sample[j]
        if(!is.na(s)&&nzchar(s))return(s)
      }
      source_column(column)
    },character(1))
    repeated<-duplicated(sample_labels)|duplicated(sample_labels,fromLast=TRUE)
    sample_labels[repeated]<-paste0(sample_labels[repeated]," [",vapply(sample_columns[repeated],source_column,character(1)),"]")
    for(column in sample_columns){
      m<-mm[mm$source_column==column,,drop=FALSE];v<-m$value[match(ids,m$peptide_id)]
      state<-m$value_status[match(ids,m$peptide_id)];v[!is.finite(v)|state %in% c("missing","invalid","non_numeric")]<-NA_real_
      label<-paste0("Quantity: ",sample_labels[match(column,sample_columns)])
      add_peptide_column(label,paste("source",source_column(column)),v)
    }
  }
  # Reserve every ordinary heading before adding contextual suffixes, so a
  # generated heading cannot take another sample's literal name. All dynamic
  # fields share this namespace, including comparison effects and quantities.
  labels<-vapply(peptide_columns,function(x)x$label,character(1))
  collisions<-labels %in% names(peptides)|duplicated(labels)|duplicated(labels,fromLast=TRUE)
  reserved<-unique(c(names(peptides),labels))
  for(i in seq_along(peptide_columns)){
    column<-peptide_columns[[i]];label<-column$label
    if(collisions[i]){
      stem<-paste0(label," [",column$context,"]");label<-stem;suffix<-2L
      while(label %in% c(reserved,names(peptides))){label<-paste0(stem," (",suffix,")");suffix<-suffix+1L}
    }
    peptides[[label]]<-column$value
  }
  other<-if(length(ids)&&nrow(ev))peptides[!peps$usage %in% "Used in eligible pair" & vapply(ids,function(id)any(ev$role[ex[[id]]] %in% c("N","C")),logical(1)),,drop=FALSE] else empty()
  annotations<-bind(list(if(is.data.frame(payload$merops$matches)&&nrow(payload$merops$matches))data.frame(Source="MEROPS",payload$merops$matches,check.names=FALSE),if(is.data.frame(payload$uniprot$site_annotations)&&nrow(payload$uniprot$site_annotations))data.frame(Source="UniProt",payload$uniprot$site_annotations,check.names=FALSE)))
  if("site_id" %in% names(annotations)){
    j<-match(annotations$site_id,sites$site_id)
    annotations$`Annotation ID`<-annotations$site_id
    annotations$Protein<-sites$accession[j]
    annotations$Bond<-paste(sites$bond_left[j],sites$bond_right[j],sep=" | ")
    leading<-c("Annotation ID","Protein","Bond","Source")
    annotations<-annotations[c(leading,setdiff(names(annotations),leading))]
  }
  settings<-list();add<-function(section,name,value){settings[[length(settings)+1L]]<<-data.frame(Section=section,Setting=name,Value=as.character(value),check.names=FALSE)}
  flatten<-function(section,name,value){
    if(is.null(value))return(invisible(NULL))
    if(is.data.frame(value)){for(i in seq_len(nrow(value)))for(n in names(value))flatten(section,paste(name,i,n,sep=" / "),value[[n]][i]);return(invisible(NULL))}
    if(is.list(value)){nn<-names(value);if(is.null(nn))nn<-as.character(seq_along(value));for(i in seq_along(value))flatten(section,paste(name,nn[i],sep=" / "),value[[i]]);return(invisible(NULL))}
    add(section,name,paste(as.character(value),collapse="; "))
  }
  if(!mapping_only)add("Method","Representative selection","Fixed digestion/coverage/overlap/intensity/ID rule; results depend on this choice; alternatives are available in All pairs and audit.")
  add("Method","Sample quantities","Values used by the report, with no additional export transformation, imputation or rounding. Wizard log2 inputs are converted to linear before analysis; sample scale columns record this conversion. Missing values are blank; numeric zero is preserved.")
  add("Export","Scope",if(mapping_only)"All supplied forms" else "All results at default display coverage")
  if(!mapping_only)add("Export","Ratio threshold",1)
  add("Export","Filters",if(mapping_only)"No filters" else "Default coverage; no search, pattern or evidence filter")
  for(n in names(payload$cleavage$profile))flatten("Digestion",n,payload$cleavage$profile[[n]])
  for(n in c("selection","inference_settings"))if(!is.null(pv[[n]]))flatten("Analysis",n,pv[[n]])
  if(is.data.frame(pv$samples))for(i in seq_len(nrow(pv$samples))){
    ss<-pv$samples[i,,drop=FALSE]
    row<-data.frame(Section="Samples",Setting="",Value="",Sample=ss$Sample,Condition=ss$Condition,
      Include=ss$Use,`Source column`=source_column(ss[["Peptide column"]]),check.names=FALSE)
    q<-payload$column_mapping_audit$quantities
    row$`Input scale`<-row$`Export scale`<-"Not recorded"
    if(is.data.frame(q)&&all(c("canonical_column","input_scale") %in% names(q))){
      j<-match(ss[["Peptide column"]],q$canonical_column)
      if(!is.na(j)){row$`Input scale`<-q$input_scale[j];row$`Export scale`<-"linear"}
    }
    if("Experimental unit" %in% names(ss))row$`Independent unit`<-ss[["Experimental unit"]]
    settings[[length(settings)+1L]]<-row
  }
  if(!is.null(payload$reference_access))flatten("References","Versions and hashes",payload$reference_access)
  if(!is.null(payload$export_settings))for(n in names(payload$export_settings))flatten("Run",n,payload$export_settings[[n]])
  if(is.data.frame(pv$comparisons))for(i in seq_len(nrow(pv$comparisons))){cp<-pv$comparisons[i,,drop=FALSE];section<-paste0("Comparison: ",cp$Name)
    for(n in names(cp))if(n!="Name")flatten(section,n,cp[[n]])
    totals<-vapply(c(cp$Control,cp$Treatment),function(cond)sum(pv$samples$Use=="Yes"&pv$samples$Condition==cond),integer(1))
    add(section,"Control total",totals[1]);add(section,"Treatment total",totals[2])
    if(!mapping_only){
      add(section,"Default display coverage k",max(1,min(3,totals)))
      ks<-pv$selection$k_by_comparison[[cp$Name]];if(is.null(ks))ks<-max(1,min(3,totals))
      add(section,"Frozen selection k",ks)
    }
  }
  overview<-display(chosen,TRUE)
  # Without a contrast there is no condition coverage threshold or effect. Keep
  # the same descriptive rows shown at k=0 in HTML, with a null comparison key.
  visible<-if("Passes default coverage" %in% names(overview))overview[overview$`Passes default coverage` | overview$Status %in% "No comparison",,drop=FALSE] else overview
  # Presentation scope only: retain all pairs and observations in audit tables.
  # One-condition joint quantification remains a result; a component absent
  # throughout the comparison does not. Do not move these rows to Other candidates.
  if("Status" %in% names(visible))
    visible<-visible[!visible$Status %in% c("Spanning quantified; neo not quantified","Neo quantified; spanning not quantified"),,drop=FALSE]
  tables<-list(Results=public_results(visible),Peptides=peptides,Settings=bind(settings),`All pairs`=display(all),`Other candidates`=other,Annotations=annotations)
  attr(tables,"overview")<-overview
  attr(tables,"pair_summaries")<-prepared
  tables
}
