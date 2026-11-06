// Report presentation contract. Arithmetic/selection behavior lives in test_pair_summary.cjs.
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const root=path.join(__dirname,'..');
const template=fs.readFileSync(path.join(root,'inst','templates','cleavage_report_template.html'),'utf8');
const helper=fs.readFileSync(path.join(root,'inst','templates','pair_summary.js'),'utf8');
const helperMarker='__ATELENTE_PAIR_SUMMARY_JS__',dataMarker='__PROTEOMASTER_DATA__';
const occurrences=(text,needle)=>text.split(needle).length-1;
const normalized=text=>text.replace(/\r\n/g,'\n').trim();
assert.equal(occurrences(template,helperMarker),1,'Template must embed the summary helper exactly once.');
assert.equal(occurrences(template,dataMarker),1,'Template must have one payload slot.');
assert.ok(!/<\/script\s*>/i.test(helper),'Embedded helper must not prematurely close its script element.');
const requiredIds=['overviewHead','sites','coverage','threshold','patternCounts','pairDetail',
  'pairChoice','pairBody','allPairBody','localTrack','peptides','referenceRows','track',
  'otherSection','otherHeading','otherBody','otherPrev','otherNext',
  'motifSection','motifDetails','motifRole','motifFamily','motifGrids','motifExport','motifCellDetail',
  'mappingChoice','mappingControls','mappingStatus'];
requiredIds.push('downloads','tableChooser','tableChoices','downloadChosen','downloadStatus');
assert.ok(!template.includes('About these results'),'Long explanation panels removed');
assert.ok(!template.includes('difference from background'),'Background correction removed');
assert.ok(template.includes('Not in supplied snapshot'),'Missing reference scope must be explicit.');
assert.ok(template.includes('excluded_statistics'),'Archived test exclusion reasons must remain available.');
assert.ok(template.includes('Ratio threshold |log2FC|'),'Threshold control must identify the ratio.');
assert.ok(!template.includes('Coverage representative: greatest shared coverage'),'Obsolete selection explanation must not be restored.');

function inspectReport(html,label){
  const scripts=[...html.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script\s*>/gi)];
  assert.ok(scripts.length>=3,label+': JSON, helper and application scripts must be self-contained.');
  const markup=html.replace(/<script\b[^>]*>[\s\S]*?<\/script\s*>/gi,'');
  const ids=[...markup.matchAll(/\bid\s*=\s*["']([^"']+)["']/g)].map(x=>x[1]);
  for(const id of requiredIds)assert.equal(ids.filter(x=>x===id).length,1,label+': expected unique '+id);
  for(const id of ['overviewMatrix','pairPanels','pairPlot','diagnosticSection','diagnosticBody'])assert.ok(!ids.includes(id),label+': obsolete view '+id);
  assert.ok(!/class=["'][^"']*\bpair-chart\b/.test(html),label+': sample charts must not replace the position map.');
  let executable=0,helpers=0,payloads=0;
  for(const [index,script] of scripts.entries()){
    assert.ok(!/\bsrc\s*=/i.test(script[1]),label+': scripts must work offline.');
    if(/\btype=["']application\/json["']/i.test(script[1])){
      assert.equal(/\bid=["']data["']/i.test(script[1]),true,label+': expected report payload.');
      require('../inst/templates/report_payload.js').parseSync(script[2]);payloads++;continue;
    }
    executable++;
    new vm.Script(script[2],{filename:label+':script-'+index});
    if(normalized(script[2])===normalized(helper)){
      helpers++;
      // Execute the actual embedded helper in the browser branch, without require/module.
      const sandbox={};vm.runInNewContext(script[2],sandbox);
      for(const name of ['summarizePairs','rankPairs','view','mixedDirections','defaultK','motifSummary'])
        assert.equal(typeof sandbox.PairSummary?.[name],'function',label+': missing helper API '+name);
    }
  }
  assert.equal(payloads,1,label+': expected one JSON payload.');
  assert.equal(helpers,1,label+': expected the complete current helper source.');
  assert.ok(executable>=2,label+': missing executable report application.');
}

// Use callback replacement so JavaScript $ sequences are preserved literally.
inspectReport(template.replace('__ATELENTE_REPORT_PAYLOAD_JS__',()=>fs.readFileSync(path.join(root,'inst','templates','report_payload.js'),'utf8')).replace('__ATELENTE_WORKBOOK_JS__',()=>fs.readFileSync(path.join(root,'inst','templates','report_workbook.js'),'utf8')).replace(helperMarker,()=>helper).replace(dataMarker,()=>JSON.stringify({})), 'template');
for(const reportPath of process.argv.slice(2))inspectReport(fs.readFileSync(reportPath,'utf8'),reportPath);
console.log('Pair-view template syntax, offline helper embedding and table/map presentation contracts passed.');
