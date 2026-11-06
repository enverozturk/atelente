'use strict';
const assert=require('node:assert/strict');
const P=require('../inst/templates/pair_summary.js');
const comparison={Name:'T-C',Control:'Control',Treatment:'Treatment'};
const samples=['c1','c2','c3','t1','t2','t3'].map(x=>({Sample:x,Condition:x[0]==='c'?'Control':'Treatment',Use:'Yes','Peptide column':x}));
samples.push({Sample:'Excluded',Condition:'Control',Use:'No','Peptide column':'excluded'},
  {Sample:'Other',Condition:'Other',Use:'Yes','Peptide column':'other'});
const preferred={neo_other_end_compatible:true,spanning_n_compatible:true,spanning_c_compatible:true,
  neo_missed_cleavages:0,spanning_missed_cleavages:0};
function rows(id,counts=[3,3],intensity=100){
  return samples.map(s=>{
    const included=s.Use==='Yes'&&['Control','Treatment'].includes(s.Condition);
    const i=Number(s.Sample.slice(1))-1,condition=s.Condition==='Control'?0:1;
    const n=included&&i<counts[condition]?intensity:null,span=included&&i<counts[condition]?intensity:null;
    return {site_id:'site',record_id:1,accession:'P1',bond_left:10,bond_right:11,role:'N',neo_peptide_id:'n'+id,spanning_peptide_id:'s'+id,
      source_column:s['Peptide column'],neo_value:n,spanning_value:span,neo_status:n===null?'missing':'observed',
      spanning_status:span===null?'missing':'observed',log2_ratio:n===null?null:0};
  });
}
function pair(id,counts=[3,3],intensity=100,digestion=preferred){
  return {...P.summarizePairs(rows(id,counts,intensity),samples,comparison)[0],digestion};
}
const selected=(pairs,options)=>P.rankPairs(pairs,options)[0];
const full=pair('full'),sparse=pair('sparse',[2,2],1e12);
assert.equal(selected([sparse,full]).key,full.key,'Eligibility precedes intensity.');
assert.equal(selected([sparse,full],{selectionK:2}).key,full.key,'Coverage precedes intensity within preferred class.');
const compatible=pair('compatible',[2,2],1),incompatible=pair('incompatible',[3,3],1e12,{...preferred,neo_other_end_compatible:false});
assert.equal(selected([incompatible,compatible],{selectionK:2}).key,compatible.key,'Preferred digestion precedes coverage among eligible pairs.');
assert.equal(selected([incompatible,compatible],{selectionK:3}).key,incompatible.key,'Digestion preference cannot override eligibility.');
const bright=pair('bright',[3,3],1000),dim=pair('dim',[3,3],100);
assert.equal(selected([dim,bright]).key,bright.key,'Intensity breaks equal-coverage ties.');
assert.equal(selected([bright,dim]).key,bright.key);
assert.equal(selected([pair('b'),pair('a')]).neo_peptide_id,'na','Stable form IDs resolve exact ties.');
assert.deepEqual(P.legacyRankPairs([dim,bright]).map(p=>p.key),[bright,dim].map(p=>p.key).sort(),'Legacy ordering ignores intensities.');
const belowA=pair('a',[1,1],1),belowZ=pair('z',[1,1],1e20);
const fallback=selected([belowZ,belowA]);
assert.equal(fallback.key,belowA.key);assert.equal(fallback.selection.fallback,true);
assert.equal(fallback.selection.eligible,false);assert.equal(fallback.selection.eligibility_reason,'below_selection_coverage');
assert.equal(P.rankPairs([full])[0].selection.selection_k,3);
for(const k of [0,1,2,3])P.view(full,k,1); // Display operations cannot alter frozen selection.
assert.equal(P.rankPairs([full])[0].selection.selection_k,3);
assert.throws(()=>P.rankPairs([full],{selectionK:1.5}));
assert.throws(()=>P.rankPairs([full],{selectionK:0}),/positive/);
assert.deepEqual(P.rankPairs([],{selectionK:0}),[]);
const noComparison=P.summarizePairs(rows('none'),samples,null)[0];
assert.equal(P.rankPairs([noComparison],{selectionK:0})[0].selection.eligible,false);
const absentCondition=P.summarizePairs(rows('absent'),samples.filter(s=>s.Condition!=='Treatment'),comparison)[0];
assert.equal(absentCondition.selection_k,1);
assert.equal(P.rankPairs([absentCondition])[0].selection.selection_k,1);
assert.throws(()=>P.rankPairs([full,{...full,selection_k:2}]));
assert.throws(()=>P.rankPairs([full,{...full,site_id:'another'}]));

const quantities=rows('score',[3,3],4);
quantities[0].neo_value=16;quantities[0].spanning_value=4;quantities[0].log2_ratio=2;
quantities[1].neo_value=1e100;quantities[1].spanning_value=null;quantities[1].spanning_status='missing';quantities[1].log2_ratio=null;
for(const r of quantities.filter(r=>['excluded','other'].includes(r.source_column))){r.neo_value=1e200;r.spanning_value=1e200;r.neo_status=r.spanning_status='observed';r.log2_ratio=0;}
const q=P.summarizePairs(quantities,samples,comparison)[0];
assert.ok(Math.abs(q.common_log2_intensity-2.2)<1e-12,'Only five common included samples contribute: (3+2+2+2+2)/5.');
const reversed=P.summarizePairs([...quantities].reverse(),[...samples].reverse(),comparison)[0];
assert.equal(q.common_log2_intensity,reversed.common_log2_intensity);
const swapped=P.summarizePairs(quantities,samples,{...comparison,Control:'Treatment',Treatment:'Control'})[0];
assert.equal(q.common_log2_intensity,swapped.common_log2_intensity);
assert.ok(Math.abs(q.log2_fc+swapped.log2_fc)<1e-12);
const asymRows=[...rows('left',[3,2],8),...rows('right',[2,3],16)];
const asym=P.summarizePairs(asymRows,samples,comparison);
const asymSwapped=P.summarizePairs(asymRows,samples,{...comparison,Control:'Treatment',Treatment:'Control'});
assert.equal(selected(asym,{selectionK:2}).key,selected(asymSwapped,{selectionK:2}).key,'Swapping condition names preserves representative identity.');
const unbalanced=samples.filter(s=>!['c2','c3'].includes(s.Sample));
const unequal=P.summarizePairs(quantities,unbalanced,comparison)[0];
assert.equal(unequal.common_log2_intensity,2.25,'Samples, not conditions, receive equal weight.');
const huge=pair('huge',[3,3],1e300),tiny=pair('tiny',[3,3],1e-300);
assert.equal(Number.isFinite(huge.common_log2_intensity),true);
assert.equal(Number.isFinite(tiny.common_log2_intensity),true);
assert.equal(selected([tiny,huge]).key,huge.key);
const invalid=rows('invalid');invalid[0].log2_ratio=999;
const invalidPair=P.summarizePairs(invalid,samples,comparison)[0];
assert.equal(invalidPair.common_log2_intensity,null);
assert.equal(P.selectionEligibility(invalidPair,1).reason,'invalid_arithmetic');
const single=rows('single');
for(const r of single.filter(r=>r.source_column.startsWith('c'))){r.neo_value=null;r.neo_status='missing';r.log2_ratio=null;}
const singlePair=P.summarizePairs(single,samples,comparison)[0];
assert.equal(singlePair.log2_fc,null);
assert.equal(selected([singlePair]).selection.eligible,true);
assert.equal(selected([singlePair]).selection.eligibility_reason,'neo_quantified_only_in_treatment');
const zero=single.map(r=>({...r}));zero[0].neo_value=0;zero[0].neo_status='zero';
assert.equal(P.selectionEligibility(P.summarizePairs(zero,samples,comparison)[0],3).eligible,false);

// The source statistical result and all effects survive selection unchanged.
const original={...bright,p_value:.4,bh_adjusted_p:.9,test_family:'original-family'};
const snapshot=JSON.stringify(original),ranked=P.rankPairs([original]);
assert.equal(JSON.stringify(original),snapshot);
for(const field of ['p_value','bh_adjusted_p','test_family','log2_fc','neo_log2_fc','spanning_log2_fc'])assert.equal(ranked[0][field],original[field]);
const effectA={...pair('effectA'),log2_fc:100,p_value:1e-100},effectB={...pair('effectB'),log2_fc:-.01,p_value:.9};
assert.equal(selected([effectB,effectA]).key,effectA.key);
assert.equal(selected([{...effectA,log2_fc:-1e9,p_value:1},{...effectB,log2_fc:1e9,p_value:0}]).key,effectA.key);

function fixture(){
  const p={site_id:'site',record_id:1,role:'N',neo_peptide_id:'n',spanning_peptide_id:'s'};
  const neo={site_id:'site',record_id:1,peptide_id:'n',role:'N',start:11,end:15,is_candidate:true};
  const span={site_id:'site',record_id:1,peptide_id:'s',role:'spanning',start:8,end:15,is_candidate:false};
  const terminal=(e,role)=>({...e,role,digestion_compatible:true,pre_digest_n_label:false,processing_status:'not_assessed'});
  return {p,payload:{peptides:[{peptide_id:'n',sequence:'AKPAR'},{peptide_id:'s',sequence:'AKPARAAA'}],
    cleavage:{profile:{enzyme:'trypsin',lysine_blocked:false},evidence:[neo,span],
      termini:[terminal(neo,'N'),terminal(neo,'C'),terminal(span,'N'),terminal(span,'C')]}}};
}
let f=fixture(),d=P.buildDigestionLookup(f.payload)(f.p);
assert.equal(d.neo_missed_cleavages,0,'K followed by P is not a missed cleavage; terminal R is excluded.');
assert.equal(d.spanning_missed_cleavages,1,'Internal R is counted.');
assert.equal(P.digestionPreference(d).preferred,true,'Internal MC is recorded but does not disqualify compatible ends.');
f.payload.peptides[1].sequence='AKPAAAAA';
assert.equal(P.digestionPreference(P.buildDigestionLookup(f.payload)(f.p)).preferred,true);
f.payload.peptides[1].sequence='AKAAAAAA';
assert.equal(P.buildDigestionLookup(f.payload)(f.p).spanning_missed_cleavages,1);
f.payload.cleavage.profile.lysine_blocked=true;
assert.equal(P.buildDigestionLookup(f.payload)(f.p).spanning_missed_cleavages,0,'Blocked K is not counted.');
f.payload.cleavage.termini[0].pre_digest_n_label=true;
assert.equal(P.buildDigestionLookup(f.payload)(f.p).neo_other_end_compatible,true,'Neo-N target N label is not its other end.');
f.payload.cleavage.termini[2].pre_digest_n_label=true;
assert.equal(P.buildDigestionLookup(f.payload)(f.p).spanning_n_compatible,false,'Pre-digest labelled spanning N is not digest-generated.');
f.payload.cleavage.termini[2].processing_status='reference_end';
assert.equal(P.buildDigestionLookup(f.payload)(f.p).spanning_n_compatible,null);
f.payload.cleavage.termini[2].processing_status='possible_initiator_met_removal';
assert.equal(P.buildDigestionLookup(f.payload)(f.p).spanning_n_compatible,null);
f.payload.cleavage.profile.enzyme='unknown';
d=P.buildDigestionLookup(f.payload)(f.p);
assert.equal(d.spanning_missed_cleavages,null);assert.equal(d.neo_other_end_compatible,null);
f=fixture();f.p.role='C';f.payload.cleavage.evidence[0].role='C';f.payload.cleavage.termini[0].pre_digest_n_label=true;
assert.equal(P.buildDigestionLookup(f.payload)(f.p).neo_other_end_compatible,false,'Neo-C other N end respects pre-digestion labelling.');
f=fixture();f.payload.cleavage.termini=f.payload.cleavage.termini.filter(t=>!(t.peptide_id==='n'&&t.role==='C'));
assert.equal(P.buildDigestionLookup(f.payload)(f.p).neo_other_end_compatible,null,'Missing exact placement annotation remains unknown.');
f=fixture();const additional={...f.payload.cleavage.evidence[1],start:5,end:12};
f.payload.cleavage.evidence.push(additional);
f.payload.cleavage.termini.push({...additional,role:'N',digestion_compatible:false,processing_status:'not_assessed'},
  {...additional,role:'C',digestion_compatible:true,processing_status:'not_assessed'});
assert.equal(P.buildDigestionLookup(f.payload)(f.p).spanning_n_compatible,false,'All possible supporting placements must meet preference.');
f.payload.cleavage.evidence.reverse();f.payload.cleavage.termini.reverse();
assert.equal(P.buildDigestionLookup(f.payload)(f.p).spanning_n_compatible,false);
assert.equal(P.buildDigestionLookup(f.payload)({...f.p,record_id:99}).neo_other_end_compatible,null,'No accession-only or other-record joins.');
assert.equal(P.digestionPreference({...preferred,spanning_missed_cleavages:null}).preferred,true);
assert.equal(P.digestionPreference({...preferred,spanning_n_compatible:null}).preferred,false);
assert.equal(P.digestionPreference({...preferred,neo_missed_cleavages:-1}).preferred,true);
for(const ambiguous of ['X','B','Z','J']){
  f=fixture();f.payload.peptides[0].sequence='AA'+ambiguous+'AA';
  assert.equal(P.buildDigestionLookup(f.payload)(f.p).neo_missed_cleavages,null);
  assert.equal(P.digestionPreference(P.buildDigestionLookup(f.payload)(f.p)).preferred,true,'Unknown MC alone does not change the independently supplied end annotations.');
}
for(const residue of ['U','O']){
  f=fixture();f.payload.peptides[0].sequence='AA'+residue+'AA';
  assert.equal(P.buildDigestionLookup(f.payload)(f.p).neo_missed_cleavages,0);
}
const lookup=new Map([[full.key,preferred]]);
assert.equal(P.rankPairs([full],{digestionByPair:lookup})[0].selection.preferred_digestion,true);
assert.equal(P.rankPairs([full],{digestionByPair:()=>({...preferred,neo_missed_cleavages:1})})[0].selection.preferred_digestion,true);
// V2 removes MC from every ranking stage, including unknown or absent counts.
assert.equal(P.SELECTION_POLICY,'ends_coverage_overlap_intensity_v3');
const withMC=pair('withMC',[3,3],1,{...preferred,neo_missed_cleavages:2,spanning_missed_cleavages:3});
const withoutMC=pair('withoutMC',[2,2],1e12);
assert.equal(selected([withoutMC,withMC],{selectionK:2}).key,withMC.key,'Coverage wins despite positive MC counts.');
const tieA=pair('a',[3,3],100),tieB=pair('b',[3,3],1000);
for(const mc of [0,1,5,null,undefined,-1,NaN]){
  const changed={...tieB,digestion:{...preferred,neo_missed_cleavages:mc,spanning_missed_cleavages:mc}};
  const actual=selected([tieA,changed]);
  assert.equal(actual.key,tieB.key,'MC cannot override intensity or stable IDs.');
  assert.deepEqual(actual.selection.digestion_reasons,[]);
  assert.equal(actual.selection.rule_version,P.SELECTION_POLICY);
}
// V3: overlap only breaks shared-coverage ties, before intensity.
const overlapPair=(id,count,coverage=[3,3],intensity=100)=>pair(id,coverage,intensity,{...preferred,other_candidate_bonds:count});
assert.equal(selected([overlapPair('many',4,[3,3],1e10),overlapPair('few',1,[3,3],1)]).neo_peptide_id,'nfew');
assert.equal(selected([overlapPair('many',4,[3,3]),overlapPair('few',0,[2,3])],{selectionK:2}).neo_peptide_id,'nmany','Minimum common coverage precedes overlap.');
assert.equal(selected([overlapPair('many',4,[3,2]),overlapPair('few',0,[2,2])],{selectionK:2}).neo_peptide_id,'nmany','Total common coverage precedes overlap.');
assert.equal(selected([overlapPair('unknown',null,[3,3],1e20),overlapPair('known',10,[3,3],1)]).neo_peptide_id,'nknown');
assert.equal(selected([overlapPair('zero',0),overlapPair('negative',-1)]).neo_peptide_id,'nzero');
assert.equal(selected([overlapPair('z',0,[1,1],1e20),overlapPair('a',9,[1,1],1)]).neo_peptide_id,'na','No-eligible fallback ignores overlap.');
f=fixture();
f.payload.cleavage.sites=[{site_id:'site',record_id:1,bond_left:10},
  {site_id:'edge_start',record_id:1,bond_left:8},{site_id:'inside',record_id:1,bond_left:12},
  {site_id:'duplicate',record_id:1,bond_left:12},{site_id:'edge_end',record_id:1,bond_left:15},
  {site_id:'other_record',record_id:2,bond_left:11}];
assert.equal(P.buildDigestionLookup(f.payload)(f.p).other_candidate_bonds,2,'Distinct same-record internal bonds only; focal and end boundary excluded.');
f.payload.cleavage.sites.reverse();
assert.equal(P.buildDigestionLookup(f.payload)(f.p).other_candidate_bonds,2);
f.payload.cleavage.sites=f.payload.cleavage.sites.filter(s=>s.site_id==='site');
assert.equal(P.buildDigestionLookup(f.payload)(f.p).other_candidate_bonds,0);
delete f.payload.cleavage.sites;
assert.equal(P.buildDigestionLookup(f.payload)(f.p).other_candidate_bonds,null,'Missing site data cannot establish zero overlap.');
console.log('Representative selection: eligibility, end compatibility, shared coverage, other-bond overlap, pooled intensity, identity and immutable source results passed.');
