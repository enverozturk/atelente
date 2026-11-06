'use strict';
const assert=require('node:assert/strict'),P=require('../inst/templates/pair_summary.js');
const positions=['P4','P3','P2','P1','P1prime','P2prime','P3prime','P4prime'];
function row(id,fc,residue='A',family='primary',role='N',co=3){
  const site={site_id:id,record_id:1,bond_left:Number(id),...Object.fromEntries(positions.map(p=>[p,residue]))};
  return {site,families:family?[family]:[],pair:{site_id:id,record_id:1,key:id,role,comparison:'T-C',
    neo_peptide_id:'neo'+id,spanning_peptide_id:'span'+id,arithmetic_valid:true,log2_fc:fc,
    control:{co,total:3},treatment:{co,total:3}}};
}
function alternative(source,record,bond){
  const clone=structuredClone(source),id=record+':'+bond;
  Object.assign(clone.site,{site_id:id,record_id:record,bond_left:bond});
  Object.assign(clone.pair,{site_id:id,record_id:record,key:id});return clone;
}
const get=(result,group,residue,position='P1')=>result.groups[group].cells.find(c=>c.residue===residue&&c.position===position);
const options={role:'N',family:'primary',k:3,tau:1};
const input=[row('1',2,'A'),row('2',-2,'S'),row('3',0,'A'),row('4',3,'C','secondary'),row('5',2,'C','primary','C'),row('6',2,'C','primary','N',1)];
input.push(input[0],alternative(input[0],2,100));
const snapshot=JSON.stringify(input),result=P.motifSummary(input,options);
assert.equal(result.measurement_pairs,3);assert.deepEqual(result.groups.map(g=>g.measurement_pairs),[1,1]);
assert.equal(get(result,0,'A').numerator,1);assert.equal(get(result,0,'A').denominator,1);
assert.equal(get(result,0,'A').percent,100);assert.equal(get(result,1,'A').percent,0);
assert.equal(get(result,1,'S').percent,100);
assert.ok(!JSON.stringify(result).match(/background|difference_pp/),'No background correction fields');
assert.deepEqual(P.motifSummary(input.slice().reverse(),options),result,'Mapping order does not select a reference context');
assert.equal(JSON.stringify(input),snapshot,'Presentation leaves source observations untouched');

// Conflicting or missing context affects only that measurement group's position.
const conflict=alternative(input[0],3,200);conflict.site.P1='G';conflict.site.P3='X';conflict.site.P4=null;
const partial=P.motifSummary([input[0],conflict,row('7',2,'C'),row('8',-2,'S')],options);
assert.equal(partial.groups[0].measurement_pairs,2);
assert.equal(get(partial,0,'A','P1').numerator,0);assert.equal(get(partial,0,'G','P1').numerator,0);
assert.equal(get(partial,0,'C','P1').denominator,1);assert.equal(get(partial,0,'C','P1').percent,100);
for(const p of ['P3','P4'])assert.equal(get(partial,0,'C',p).denominator,1);
assert.equal(get(partial,0,'A','P2').denominator,2);assert.equal(get(partial,0,'A','P2').percent,50);
assert.equal(get(partial,0,'C','P2').percent,50);
const allMissing=P.motifSummary([input[0],conflict],options);
assert.equal(get(allMissing,0,'A','P1').denominator,0);assert.equal(get(allMissing,0,'A','P1').percent,null);
assert.deepEqual(P.motifSummary([input[0],conflict,row('7',2,'C'),row('8',-2,'S')].reverse(),options),partial);

// Identical percentages mean identical colors even when panel sample counts differ.
const unequal=P.motifSummary([row('1',2,'A'),row('2',2,'A'),row('3',2,'S'),row('4',2,'S'),row('5',-2,'A'),row('6',-2,'S')],options);
assert.deepEqual(unequal.groups.map(g=>g.measurement_pairs),[4,2]);
assert.equal(get(unequal,0,'A').numerator,2);assert.equal(get(unequal,1,'A').numerator,1);
assert.equal(get(unequal,0,'A').percent,50);assert.equal(get(unequal,1,'A').percent,50);

// Form IDs (not sequences), role, comparison and selected family bound the unit.
const otherForm=structuredClone(input[0]);otherForm.pair.neo_peptide_id='different_modification_form';
assert.equal(P.motifSummary([input[0],otherForm],options).groups[0].measurement_pairs,2);
const otherComparison=structuredClone(input[0]);otherComparison.pair.comparison='T2-C';
assert.equal(P.motifSummary([input[0],otherComparison],options).groups[0].measurement_pairs,2);
assert.equal(P.motifSummary(input,{...options,family:'secondary'}).groups[0].measurement_pairs,1);
assert.equal(P.motifSummary(input,{...options,role:'C'}).groups[0].measurement_pairs,1);
const unrelated=alternative(input[0],4,300);unrelated.families=['secondary'];unrelated.site.P1='Y';
assert.equal(get(P.motifSummary([input[0],unrelated],options),0,'A').percent,100,'Other families cannot invalidate context');
const invalid=row('9',2);invalid.pair.arithmetic_valid=false;
assert.equal(P.motifSummary([invalid],options).measurement_pairs,0);

const zero=P.motifSummary([row('1',0),row('2',0.1),row('3',-0.1)],{...options,tau:0});
assert.deepEqual(zero.groups.map(g=>g.measurement_pairs),[1,1]);assert.equal(zero.measurement_pairs,3);
const empty=P.motifSummary(input,{...options,family:'absent'});
assert.equal(empty.measurement_pairs,0);assert.ok(empty.groups.every(g=>g.cells.every(c=>c.percent===null)));
assert.equal(P.motifSummary([row('1',1,'A','')],{role:'N',family:''}).measurement_pairs,1);
assert.throws(()=>P.motifSummary([row('1',1),row('1',-1)],options),/Conflicting motif quantities/);
const badCoverage=alternative(input[0],2,100);badCoverage.pair.control.co=4;badCoverage.pair.control.total=4;
assert.throws(()=>P.motifSummary([input[0],badCoverage],options),/Conflicting motif quantities/);
const badContext=structuredClone(input[0]);badContext.site.P1='Y';
assert.throws(()=>P.motifSummary([input[0],badContext],options),/Conflicting motif contexts/);
const missingID=structuredClone(input[0]);delete missingID.pair.neo_peptide_id;
assert.throws(()=>P.motifSummary([missingID],options),/peptide-form identifiers/);
assert.throws(()=>P.motifSummary(input,{role:'both'}),/one terminal role/);
assert.throws(()=>P.motifSummary(input,{role:'N',tau:-1}),/Invalid/);
console.log('Motif summary: measurement-pair deduplication, mapping consensus, family/role/coverage scope, positional missingness and observed percentages passed.');
