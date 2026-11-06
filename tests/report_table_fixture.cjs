// Explicitly synthetic presentation fixture; not an experimental validation dataset.
// node tests/report_table_fixture.cjs NEW_OUTPUT_DIRECTORY
'use strict';
const fs=require('node:fs'),path=require('node:path');
function fixture(){
  const samples=['C1','C2','C3','C4','C5','T1','T2','T3','T4','T5'].map(s=>({Sample:s,Condition:s[0]==='C'?'Control':'Treatment',Use:'Yes','Peptide column':s,'Protein column':''}));
  const D={scope:'reported',fasta_records:8,dataset_label:'Synthetic interface fixture · not experimental results',
    summary:[],proteins:[],peptides:[],matches:[],
    cleavage:{sites:[],evidence:[],termini:[],profile:{enzyme:'trypsin',lysine_blocked:false,pre_digest_n_label:''}},
    pair_view:{samples,comparisons:[{Name:'T_vs_C',Control:'Control',Treatment:'Treatment'}],ratios:[],contrasts:[],measurements:[],statistics:[]},
    merops:{matches:[]},uniprot:{records:[],site_annotations:[]}};
  function addSite(id,label){
    const accession='FIXTURE_'+id,site_id='synthetic:'+id+':20';
    D.proteins.push({record_id:id,accession,header:label,length:80,gene_symbol:label});
    D.uniprot.records.push({record_id:id,status:'matched'});
    const site={site_id,record_id:id,accession,bond_left:20,bond_right:21,context_display:'AAAA|GGGG',context_complete:true,neo_n_count:0,neo_c_count:0,spanning_count:0};
    D.cleavage.sites.push(site);return site;
  }
  function addPeptide(site,id,role,modifications=''){
    const start=role==='N'?21:role==='C'?13:15,end=role==='N'?28:role==='C'?20:25;
    D.peptides.push({peptide_id:id,source_row:D.peptides.length+2,sequence:role==='N'?'GGGGAAAA':role==='C'?'AAAAGGGG':'AAAAAAGGGGG',modifications,modifications_available:true,mapping_status:'mapped',match_count:1,exclusion_reason:'',contaminant_screening:''});
    const placement={placement_id:'placement_'+id,peptide_id:id,record_id:site.record_id,accession:site.accession,start,end,protein_length:80,is_master:false,reported_by_pd:false,pd_coordinate_agrees:null};
    D.matches.push(placement);
    const evidence={...placement,site_id:site.site_id,role,bond_left:20,bond_right:21,internal_bond:true,digestion_compatible:role==='spanning',pre_digest_n_label:'',processing_status:'not_assessed',classification:role==='spanning'?'spanning':'neo_'+role+'_candidate',is_candidate:role!=='spanning',candidate_reason:'fixture'};
    D.cleavage.evidence.push(evidence);
    for(const endRole of ['N','C'])D.cleavage.termini.push({...evidence,role:endRole,classification:endRole==='N'?'digestion_compatible':'non_digestion_compatible',digestion_compatible:endRole==='N'});
    if(role==='spanning')site.spanning_count++;else site[role==='N'?'neo_n_count':'neo_c_count']++;
    return id;
  }
  function addPair(site,id,control,treatment,mod=''){
    const neo=addPeptide(site,id+'_neo','N',mod),span=addPeptide(site,id+'_span','spanning');
    for(let i=0;i<samples.length;i++){
      const [n,s]=(i<5?control:treatment)(i%5),sample=samples[i];
      const status=v=>v===null?'missing':v===0?'zero':'observed';
      D.pair_view.ratios.push({site_id:site.site_id,record_id:site.record_id,accession:site.accession,bond_left:20,bond_right:21,role:'N',neo_peptide_id:neo,spanning_peptide_id:span,source_column:sample.Sample,sample:sample.Sample,abundance_type:'fixture',neo_value:n,spanning_value:s,neo_status:status(n),spanning_status:status(s),ratio:n>0&&s>0?n/s:null,log2_ratio:n>0&&s>0?Math.log2(n)-Math.log2(s):null});
      for(const [peptide_id,value] of [[neo,n],[span,s]])D.pair_view.measurements.push({peptide_id,source_column:sample.Sample,value,value_status:status(value),included_in_analysis:true});
    }
  }
  const constant=(n,s)=>()=>[n,s];
  const mixed=addSite(1,'MIXED');
  addPair(mixed,'mixed_a',constant(10,10),constant(40,10),'Literal <tag> & "quotes"');
  addPair(mixed,'mixed_b',constant(40,10),constant(10,10),'1xOxidation [M1]');
  addPair(addSite(2,'SPARSE_1_OF_5'),'sparse',i=>i===0?[10,10]:[null,null],i=>i===0?[40,10]:[null,null]);
  addPair(addSite(3,'NEO_ONLY_TREATMENT'),'neo_only',constant(null,10),constant(40,10));
  addPair(addSite(4,'BOTH_ONLY_TREATMENT'),'both_only',constant(null,null),constant(40,10));
  addPair(addSite(5,'SPANNING_ONLY'),'span_only',constant(null,10),constant(null,20));
  const noSpan=addSite(6,'NO_SPANNING');addPeptide(noSpan,'no_span_neo','C');
  for(const s of samples)D.pair_view.measurements.push({peptide_id:'no_span_neo',source_column:s.Sample,value:20,value_status:'observed',included_in_analysis:true});
  addPair(addSite(7,'NO_QUANTITIES'),'no_quantity',constant(null,null),constant(null,null));
  addPair(addSite(8,'SMALL_OPPOSITE_NOT_MIXED'),'small_a',constant(10,10),constant(40,10));
  addPair(D.cleavage.sites[7],'small_b',constant(10,10),constant(9.8,10));
  D.merops.matches=[{site_id:mixed.site_id,relation_type:'reference_site_candidate',context_status:'agreement',substrate_accession:mixed.accession,substrate_name:'Fixture substrate',p1_position:20,p4_p4prime:'AAAAGGGG',protease_name:'Fixture protease exact',protease_id:'FIX01',merops_record_id:'SYN_EXACT',reference_ids:'Fixture only'},
    {site_id:mixed.site_id,relation_type:'motif_compatibility',context_status:'agreement',substrate_accession:'OTHER_FIXTURE',substrate_name:'Other fixture',p1_position:77,p4_p4prime:'AAAAGGGG',protease_name:'Fixture protease motif',protease_id:'FIX02',merops_record_id:'SYN_MOTIF',reference_ids:'Fixture only'}];
  D.uniprot.site_annotations=[{site_id:mixed.site_id,record_id:1,feature_type:'CHAIN',boundary:'start',bond_left:20,bond_right:21,begin:21,end:80,description:'Synthetic annotation <literal>',evidence_codes:'Fixture only'}];
  return D;
}
function writeReport(D,file){
  const root=path.join(__dirname,'..');
  let html=fs.readFileSync(path.join(root,'inst','templates','cleavage_report_template.html'),'utf8');
  html=html.replace('__ATELENTE_REPORT_PAYLOAD_JS__',()=>fs.readFileSync(path.join(root,'inst','templates','report_payload.js'),'utf8'));
  html=html.replace('__ATELENTE_WORKBOOK_JS__',()=>fs.readFileSync(path.join(root,'inst','templates','report_workbook.js'),'utf8'));
  html=html.replace('__ATELENTE_PAIR_SUMMARY_JS__',()=>fs.readFileSync(path.join(root,'inst','templates','pair_summary.js'),'utf8'));
  html=html.replace('__PROTEOMASTER_DATA__',()=>JSON.stringify(D).replaceAll('<','\\u003c').replaceAll('&','\\u0026'));
  fs.writeFileSync(file,html,'utf8');
}
if(require.main===module){
  if(process.argv.length!==3)throw Error('NEW_OUTPUT_DIRECTORY required');
  const out=process.argv[2];fs.mkdirSync(out,{recursive:true});const D=fixture();
  writeReport(D,path.join(out,'synthetic_table_report.html'));
  delete D.pair_view;delete D.merops;delete D.uniprot;D.dataset_label+=' · legacy mapping-only payload';
  writeReport(D,path.join(out,'legacy_mapping_only.html'));
  console.log('Synthetic and legacy UI fixtures generated in '+path.resolve(out));
}
module.exports={fixture,writeReport};
