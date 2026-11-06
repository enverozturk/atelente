'use strict';
const assert=require('node:assert/strict'),P=require('../inst/templates/pair_summary.js');
function pair(id,fc,site='s',co=3){return {site_id:site,record_id:1,key:site+id,comparison:'T-C',role:'N',neo_peptide_id:id,spanning_peptide_id:'span',log2_fc:fc,neo_log2_fc:fc,spanning_log2_fc:0,common_log2_intensity:3,arithmetic_valid:true,control:{co,total:3,neo:co,spanning:co},treatment:{co,total:3,neo:co,spanning:co}};}
const ps=[...Array.from({length:8},(_,i)=>pair('up'+i,2)),pair('down1',-2),pair('down2',-2),pair('small',-.01),pair('sparse',-3,'s',1)];
assert.equal(P.directionSummary(ps,3,1).label,'Predominantly increased');
assert.equal(P.directionSummary(ps,3,1).below_threshold,1);
assert.equal(P.directionSummary(ps.slice(1),3,1).label,'Mixed directions');
assert.equal(P.directionSummary(ps.slice(0,8),3,1).label,'Increased');
assert.equal(P.directionSummary([pair('d',-2),pair('small',.01)],3,1).label,'Decreased');
assert.equal(P.directionSummary([pair('a',.01)],3,1).label,'Below threshold');
assert.equal(P.directionSummary([pair('a',0)],3,0).directional,0);
assert.equal(P.directionSummary([pair('a',.01),pair('b',-.01)],3,0).label,'Mixed directions');
assert.equal(P.directionSummary([pair('a',1,'s',1)],3,1).label,'No computable pairs');
assert.equal(P.directionSummary([pair('a',2),pair('a',2,'mapping2')],3,1).increased,1,'Alternative mapping does not multiply directions');
assert.throws(()=>P.directionSummary([pair('a',2),pair('a',-2,'mapping2')],3,1),/Conflicting/);
function row(p){return {key:p.site_id,site:{site_id:p.site_id,record_id:p.record_id,bond_left:20},role:p.role,selected:p,pairs:[p]};}
const a=pair('a',2),b={...pair('a',2,'alternate'),record_id:2};
const input=[row(a),row(b),row(pair('different_form',2))],snapshot=JSON.stringify(input),groups=P.groupOverviewRows(input);
assert.equal(groups.length,2);assert.equal(groups[0].mappings.length,2);assert.equal(groups[0].pairs.length,1);assert.equal(groups[0].computable,1);
assert.equal(JSON.stringify(input),snapshot,'Grouping leaves source rows unchanged');
assert.equal(P.groupOverviewRows([row(a),row({...b,role:'C'})]).length,2);
assert.equal(P.groupOverviewRows([row(a),row({...b,comparison:'Other'})]).length,2);
assert.equal(P.groupOverviewRows([row(a),row({...b,neo_peptide_id:'span',spanning_peptide_id:'a'})]).length,2,'Directed pairs remain distinct');
assert.throws(()=>P.groupOverviewRows([row(a),row({...b,log2_fc:3})]),/inconsistent/);
console.log('Overview grouping and direction summaries: directed IDs, mappings retained, measurement deduplication, 80% boundary, zero/coverage rules passed.');
