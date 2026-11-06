const assert=require('node:assert/strict');
const P=require('../inst/templates/pair_summary.js');
const samples=['c1','c2','c3','t1','t2','t3'].map(x=>({Sample:x,Condition:x[0]==='c'?'Control':'Treatment',Use:'Yes','Peptide column':x}));
samples.push({Sample:'Excluded',Condition:'Control',Use:'No','Peptide column':'x'});
const comparison={Name:'T-C',Control:'Control',Treatment:'Treatment'};
function row(sample,n,s,neo='n',span='s',role='N',site='site'){
  const status=x=>x===null?'missing':x===0?'zero':x<0?'negative':'observed';
  return {site_id:site,role,record_id:1,accession:'P1',bond_left:10,bond_right:11,neo_peptide_id:neo,spanning_peptide_id:span,
    source_column:sample,sample,neo_value:n,spanning_value:s,neo_status:status(n),spanning_status:status(s),
    log2_ratio:n>0&&s>0?Math.log2(n)-Math.log2(s):null};
}
const rows=[row('c1',10,10),row('c2',null,100),row('c3',40,null),row('t1',80,20),row('t2',80,20),row('t3',null,null),row('x',1e9,1)];
const p=P.summarizePairs(rows,samples,comparison)[0];
assert.equal(p.control.co,1);assert.equal(p.control.total,3);assert.equal(p.treatment.co,2);
assert.ok(Math.abs(p.log2_fc-2)<1e-12);assert.ok(Math.abs(p.neo_log2_fc-3)<1e-12);assert.ok(Math.abs(p.spanning_log2_fc-1)<1e-12);
assert.ok(Math.abs(p.neo_log2_fc-p.spanning_log2_fc-p.log2_fc)<1e-12);
assert.notEqual(p.neo_all_log2_fc-p.spanning_all_log2_fc,p.log2_fc);
assert.equal(P.view(p,1,1).coverage_pass,true);assert.equal(P.view(p,3,1).coverage_pass,false);
assert.equal(P.view(p,0,1).coverage_pass,true);assert.equal(p.log2_fc,2);
assert.equal(P.defaultK(samples,comparison),3);
assert.equal(P.defaultK(samples.slice(0,1).concat(samples.slice(3,4)),comparison),1);
assert.equal(P.defaultK(samples,null),0);
assert.deepEqual(P.summarizePairs([...rows].reverse(),samples,comparison)[0].control,p.control);
const make=(control,treatment)=>['c1','c2','c3','t1','t2','t3'].map((s,i)=>row(s,...(i<3?control[i]:treatment[i-3])));
const only=P.summarizePairs(make([[null,10],[null,10],[null,null]],[[20,10],[40,10],[null,null]]),samples,comparison)[0];
assert.equal(only.log2_fc,null);assert.equal(P.view(only,2).status,'neo_quantified_only_in_treatment');
assert.equal(P.view(only,2).coverage_pass,true);assert.equal(P.view(only,3).coverage_pass,false);
const both=P.summarizePairs(make([[null,null],[null,null],[null,null]],[[20,10],[40,10],[null,null]]),samples,comparison)[0];
assert.equal(P.view(both,2).status,'both_quantified_only_in_treatment');assert.equal(P.view(both,2).coverage_pass,true);
const span=P.summarizePairs(make([[null,10],[null,10],[null,null]],[[null,20],[null,40],[null,null]]),samples,comparison)[0];
assert.equal(P.view(span,2).status,'spanning_only');assert.equal(P.view(span,2).coverage_pass,true);
const zero=P.summarizePairs(make([[0,10],[null,10],[null,null]],[[20,10],[40,10],[null,null]]),samples,comparison)[0];
assert.equal(P.view(zero,0).status,'insufficient_data');assert.equal(zero.control.neoZero,1);
const neg=P.summarizePairs(make([[-2,10],[null,10],[null,null]],[[20,10],[40,10],[null,null]]),samples,comparison)[0];
assert.equal(P.view(neg,0).status,'insufficient_data');
const rawNeo=P.summarizePairs(make([[null,10],[null,10],[10,null]],[[20,10],[40,10],[null,null]]),samples,comparison)[0];
assert.equal(P.view(rawNeo,0).status,'insufficient_data'); // neo outside co-quantified samples prevents Only in X
const copy=(base,id,effect)=>({...base,key:id,neo_peptide_id:id,log2_fc:effect});
assert.equal(P.mixedDirections([copy(p,'a',1.5),copy(p,'b',-.03)],1,1),false);
assert.equal(P.mixedDirections([copy(p,'a',1.5),copy(p,'b',-1)],1,1),true);
assert.equal(P.mixedDirections([copy(p,'a',0),copy(p,'b',-1)],1,0),false);
assert.equal(P.mixedDirections([copy(p,'a',.01),copy(p,'b',-.01)],1,0),true);
assert.equal(P.mixedDirections([copy(p,'a',1.5),copy(p,'b',-1)],3,1),false);
assert.throws(()=>P.mixedDirections([p,{...p,role:'C'}],0,1));
const tied=[copy(p,'b',999),copy(p,'a',-999)];
assert.equal(P.rankPairs(tied)[0].neo_peptide_id,'a');
assert.equal(P.componentPattern(0,-2,1),'neo_below_threshold_spanning_down');
assert.equal(P.componentPattern(-2,0,1),'neo_down_spanning_below_threshold');
assert.equal(P.componentPattern(0,0,0),'below_threshold');
// Main descriptive labels use the ratio itself, not separate component cutoffs.
const effectPair=(neo,span)=>({...p,neo_log2_fc:neo,spanning_log2_fc:span,log2_fc:neo-span});
const labels=[
  [effectPair(.8,-.8),1,'ratio_increased','Ratio increased'],
  [effectPair(-.8,.8),1,'ratio_decreased','Ratio decreased'],
  [effectPair(2.5,1),1,'ratio_increased','Ratio increased'],
  [effectPair(-2.5,-1),1,'ratio_decreased','Ratio decreased'],
  [effectPair(2,2),1,'below_ratio_threshold','Below ratio threshold'],
  [effectPair(1,0),1,'ratio_increased','Ratio increased'],
  [effectPair(-1,0),1,'ratio_decreased','Ratio decreased'],
  [effectPair(0,0),0,'below_ratio_threshold','Below ratio threshold'],
  [effectPair(.001,0),0,'ratio_increased','Ratio increased'],
  [effectPair(-.001,0),0,'ratio_decreased','Ratio decreased']
];
for(const [pair,tau,pattern,label] of labels){
  const before=JSON.stringify(pair),v=P.view(pair,1,tau);
  assert.equal(v.pattern,pattern);assert.equal(v.pattern_label,label);
  assert.equal(v.status,'quantified_in_both');assert.equal(v.coverage_pass,true);
  assert.equal(JSON.stringify(pair),before,'Display classification must not mutate pair values');
  assert.equal(P.view(pair,3,tau).coverage_pass,false,'Ratio threshold must not relax coverage');
}
const thresholdPair=effectPair(.8,-.8),originalEffect=JSON.stringify(thresholdPair);
assert.equal(P.view(thresholdPair,1,1).pattern,'ratio_increased');
assert.equal(P.view(thresholdPair,1,2).pattern,'below_ratio_threshold');
assert.equal(JSON.stringify(thresholdPair),originalEffect);
for(const unchanged of [only,both,span,zero,neg,rawNeo]){
  const a=P.view(unchanged,2,0),b=P.view(unchanged,2,10);
  assert.deepEqual(a,b,'Threshold cannot reclassify missing/one-condition states');
}
const exported={...rows[0],comparison:'T-C',log2_fc:2,fold_change:4};
const authoritative=P.summarizePairs(rows,samples,comparison,[exported])[0];
assert.equal(authoritative.fold_change,4);assert.equal(authoritative.arithmetic_valid,true);
const corrupt=P.summarizePairs(rows,samples,comparison,[{...exported,log2_fc:3}])[0];
assert.equal(corrupt.arithmetic_valid,false);assert.equal(corrupt.log2_fc,null);
assert.equal(P.view(corrupt,0).status_label,'Inconsistent source values');
assert.throws(()=>P.summarizePairs([...rows,rows[0]],samples,comparison));
assert.equal(P.view(P.summarizePairs(rows,samples,null)[0],0).status,'no_comparison');
console.log('Pair summary: shared subsets, eligibility, missingness, ranking, threshold patterns and export validation passed.');
// Optional local acceptance reports; no private data is required by the suite.
for(const filename of process.argv.slice(2)){
  const html=require('node:fs').readFileSync(filename,'utf8');
  const d=require('../inst/templates/report_payload.js').parseSync(html.match(/<script type="application\/json" id="data">([\s\S]*?)<\/script>/)[1]);
  let checked=0;
  for(const comparison of d.pair_view.comparisons){
    const summaries=P.summarizePairs(d.pair_view.ratios,d.pair_view.samples,comparison,d.pair_view.contrasts);
    for(const summary of summaries){assert.equal(summary.arithmetic_valid,true,JSON.stringify(summary.diagnostics));checked++;}
  }
  console.log('Verified '+checked+' exported pair summaries: '+filename);
}
