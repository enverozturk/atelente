// Execute application code against a small in-memory DOM. This is not browser/layout QA.
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const {fixture}=require('./report_table_fixture.cjs');
const template=fs.readFileSync(path.join(__dirname,'../inst/templates/cleavage_report_template.html'),'utf8');
class Element {
  constructor(tag){this.tagName=tag;this.children=[];this.dataset={};this.style={};this.attributes={};this.events={};this._text='';this._value=undefined;this.hidden=false;this.classList={toggle(){}};}
  append(...items){this.children.push(...items);}
  replaceChildren(...items){this.children=[...items];this._text='';if(this.tagName==='select')this._value=undefined;}
  set textContent(value){this._text=String(value);this.children=[];}
  get textContent(){return this._text+this.children.map(c=>c.textContent).join('');}
  set value(value){this._value=String(value);}
  get value(){return this._value??(this.tagName==='select'?this.children[0]?.value:'')??'';}
  setAttribute(key,value){this.attributes[key]=String(value);if(key.startsWith('data-'))this.dataset[key.slice(5)]=String(value);}
  addEventListener(name,handler){this.events[name]=handler;}
  scrollIntoView(){}
  click(){this.onclick?.();}
}
async function run(D,html=template){
  const markup=html.replace(/<script\b[^>]*>[\s\S]*?<\/script>/g,'');
  const ids=new Map();for(const match of markup.matchAll(/<([a-z][a-z0-9]*)\b[^>]*\bid="([^"]+)"[^>]*>/g)){
    const e=new Element(match[1]);const value=match[0].match(/\bvalue="([^"]*)"/);if(value)e.value=value[1];ids.set(match[2],e);
  }
  ids.set('data',new Element('script'));
  for(const match of markup.matchAll(/<select\b[^>]*id="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g))for(const option of match[2].matchAll(/<option\b[^>]*value="([^"]*)"[^>]*>([^<]*)<\/option>/g)){
    const e=new Element('option');e.value=option[1];e.textContent=option[2];ids.get(match[1]).append(e);
  }
  const errors=[];const document={getElementById(id){assert.ok(ids.has(id),'Unknown element '+id);return ids.get(id);},createElement:t=>new Element(t),createElementNS:(_,t)=>new Element(t),querySelectorAll:()=>[]};
  const scripts=[...html.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/g)],app=scripts.at(-1)[1];
  const context={document,PairSummary:require('../inst/templates/pair_summary.js'),ReportWorkbook:require('../inst/templates/report_workbook.js'),ReportPayload:{parse:async()=>D},console:{error:e=>errors.push(e)},URL,Blob,setTimeout,TextEncoder,TextDecoder};
  if(html!==template){vm.runInNewContext(scripts.at(-2)[1],context,{filename:'embedded-pair-summary.js'});vm.runInNewContext(scripts.at(-3)[1],context,{filename:'embedded-workbook.js'});assert.equal(typeof context.PairSummary?.motifSummary,'function');assert.equal(typeof context.ReportWorkbook?.mount,'function');}
  const workbook=context.ReportWorkbook;
  context.ReportWorkbook={...workbook,mount(...args){ids.exportView=args[2]?.getView;return workbook.mount(...args);}};
  await vm.runInNewContext(app,context,{filename:'report-app.js'});
  assert.deepEqual(errors,[],'Report application startup errors');assert.equal(ids.get('notice').hidden,true);
  return ids;
}
async function inspectFile(file,currentTemplate=false){
  const html=fs.readFileSync(file,'utf8'),payload=html.match(/<script\b[^>]*type="application\/json"[^>]*>([\s\S]*?)<\/script>/);
  assert.ok(payload,'Embedded JSON payload in '+file);
  const D=require('../inst/templates/report_payload.js').parseSync(payload[1]),ids=await run(D,currentTemplate?template:html);
  assert.equal(ids.get('motifGrids').children.length,2,'Two motif panels');
  for(const role of ['N','C']){ids.get('motifRole').value=role;ids.get('motifRole').onchange();assert.equal(ids.get('motifGrids').children.length,2);}
  const row=ids.get('sites').children.find(r=>typeof r.onclick==='function')??ids.get('otherBody').children.find(r=>typeof r.onclick==='function');
  if(row){row.onclick();assert.ok(ids.get('localTrack').children.length,'Selected position map');assert.ok(ids.get('peptides').children.length,'Supporting rows');
    const options=ids.get('pairChoice').children;
    if(options.length){assert.ok(options.every(o=>o.title.includes('input row')));assert.ok(options.every(o=>!o.textContent.includes('row ')));}
    if(options.length>1){ids.get('pairChoice').value=options[1].value;ids.get('pairChoice').onchange();assert.ok(ids.get('pairBody').children.length,'Alternative selected');}
  }
  const links=ids.get('downloads').children.filter(a=>a.tagName==='a'),hrefs=links.map(a=>a.href);assert.equal(new Set(hrefs).size,hrefs.length,'Unique downloads');
  for(const value of [D.excel_file,D.audit_file].filter(Boolean))assert.ok(hrefs.includes(value),'Declared primary download');
  assert.ok(!hrefs.some(h=>h.endsWith('.csv')),'Raw audit files do not clutter download UI');
  for(const href of hrefs){assert.ok(!/^[a-z]+:/i.test(href),'Local report download');assert.ok(fs.existsSync(path.resolve(path.dirname(file),href)),'Existing download '+href);}
  if(D.export_workbook?.internal_sheets?.includes('__Results')){
    assert.equal(ids.get('exportSites').textContent,'Download current view — Excel');
    const chip=ids.get('patternCounts').children.find(e=>e.textContent.startsWith('Ratio increased'));
    if(chip)chip.onclick();
    const view=ids.exportView();assert.ok(view.rows.every(r=>!chip||r.pattern==='Ratio increased'));
    assert.equal(view.rows.length,Number(ids.get('count').textContent.match(/^\d+/)[0]),'Export contains all filtered rows, not just one page');
    const dir=fs.mkdtempSync(path.join(require('node:os').tmpdir(),'atelente-view-test-'));
    fs.writeFileSync(path.join(dir,'filtered.xlsx'),await require('../inst/templates/report_workbook.js').select(D.export_workbook,['Results','Peptides','Settings'],{view}));
    fs.writeFileSync(path.join(dir,'expected.json'),JSON.stringify(view));
    console.log(JSON.stringify({filtered_export:dir,rows:view.rows.length,pattern:view.settings.pattern,threshold:view.settings.threshold,coverage:view.settings.coverage}));
  }
  console.log(JSON.stringify({report:path.resolve(file),app:currentTemplate?'current source template':'embedded report application',sites:D.cleavage.sites.length,overview:ids.get('count').textContent,download_links:links.length,detail_selected:Boolean(row),status:'in-memory DOM passed; no browser/layout test'}));
}
(async()=>{
  const D=fixture();const positions=['P4','P3','P2','P1','P1prime','P2prime','P3prime','P4prime'];
  for(const site of D.cleavage.sites)positions.forEach((p,i)=>site[p]=i<4?'A':'G');
  const ids=await run(D);
  for(const state of [ids,await run((()=>{
    const d=structuredClone(D);
    for(const r of d.pair_view.ratios.filter(r=>r.neo_peptide_id==='span_only_neo')){
      [r.neo_value,r.spanning_value]=[r.spanning_value,r.neo_value];
      [r.neo_status,r.spanning_status]=[r.spanning_status,r.neo_status];
    }
    return d;
  })())]){
    for(const k of ['3','0']){
      state.get('coverage').value=k;state.get('coverage').oninput();
      assert.ok(!state.get('sites').textContent.includes('SPANNING_ONLY'),'Missing component excluded even without coverage filter');
      assert.ok(!state.get('otherBody').textContent.includes('SPANNING_ONLY'),'Not moved to Other candidates');
      assert.ok(!state.get('patternCounts').textContent.includes('· neo not quantified'));
      assert.ok(!state.get('patternCounts').textContent.includes('· spanning not quantified'));
      assert.ok(state.get('patternCounts').textContent.includes('Quantified in one condition'),'One-condition results retained');
      assert.ok(state.exportView().rows.every(r=>!['spanning_only','neo_only'].includes(r.pattern)));
    }
  }
  assert.equal(ids.get('motifGrids').children.length,2);
  const descendants=e=>[e,...e.children.flatMap(descendants)],matrixNodes=descendants(ids.get('motifGrids'));
  assert.ok(matrixNodes.some(e=>e.style.backgroundColor==='rgb(8,81,156)'),'High-frequency cells use the requested blue');
  assert.ok(matrixNodes.filter(e=>e.tagName==='button'&&e.title?.includes('valid pairs')).every(e=>e.attributes['aria-label']===e.title),'Heatmap count/denominator/percent remains accessible');
  assert.equal(matrixNodes.filter(e=>e.className==='motif-colorbar').length,2,'Both panels use a visible shared-scale legend');
  const varied=structuredClone(D);
  varied.cleavage.sites.forEach((site,i)=>positions.forEach(p=>site[p]='ACDEFGHI'[i]));
  const variedIds=await run(varied),variedNodes=descendants(variedIds.get('motifGrids'));
  const legends=variedNodes.filter(e=>e.className==='motif-colorbar');
  assert.equal(legends[0].attributes['aria-label'],legends[1].attributes['aria-label'],'Panels share one automatic maximum');
  const upper=Number(legends[0].attributes['aria-label'].match(/to (.*) percent/)[1]);
  assert.ok(upper>0&&upper<100,'Fixture maximum is below the old fixed 100 percent scale');
  assert.ok(variedNodes.some(e=>e.style.backgroundColor==='rgb(8,81,156)'),'Observed maximum reaches darkest blue without changing cell percentages');
  assert.ok(!ids.has('motifCaption'),'Long motif explanation removed');
  assert.equal(ids.get('motifExport').disabled,false);
  assert.ok(!ids.get('downloads').children.some(a=>a.href==='neo_spanning_ratios.csv'),'Embedded measurements alone do not establish an exported file');
  assert.equal(ids.get('downloads').children.find(e=>e.tagName==='button').disabled,true,'Missing workbook disables custom download');
  const downloadable=fixture();downloadable.excel_file='atelente_results.xlsx';downloadable.audit_file='atelente_audit.zip';
  downloadable.export_workbook={base64:'fixture',sheets:['Results','Peptides','Settings','All pairs','Other candidates','Annotations']};
  const downloadIds=await run(downloadable),actions=downloadIds.get('downloads').children;
  assert.equal(actions.length,3,'Only three download actions');assert.equal(actions[0].href,'atelente_results.xlsx');assert.equal(actions[2].href,'atelente_audit.zip');
  assert.equal(actions[1].disabled,false);actions[1].onclick();assert.equal(downloadIds.get('tableChooser').hidden,false);
  const checks=downloadIds.get('tableChoices').children.map(label=>label.children[0]);
  assert.deepEqual(checks.filter(c=>c.checked).map(c=>c.value),['Results','Peptides','Settings']);
  for(const check of checks){check.checked=false;check.onchange();}assert.equal(downloadIds.get('downloadChosen').disabled,true,'Empty sheet selection disabled');
  checks[3].checked=true;checks[3].onchange();assert.equal(downloadIds.get('downloadChosen').disabled,false,'Optional table can be selected alone');
  const motif=ids.get('motifGrids').textContent;
  ids.get('search').value='DOES_NOT_EXIST';ids.get('search').events.input();
  assert.equal(ids.get('motifGrids').textContent,motif,'Search does not silently change the motif subset');
  ids.get('search').value='';ids.get('search').events.input();
  ids.get('sites').children[0].onclick();
  assert.ok(ids.get('pairChoice').children.length);
  assert.ok(ids.get('pairChoice').children.every(c=>!c.textContent.includes('row ')),'Input rows must not clutter pair options');
  assert.ok(ids.get('pairChoice').children.every(c=>c.title.includes('input row')),'Pair options retain audit identity');
  assert.ok(ids.get('localTrack').children.length,'Position map remains present');
  assert.ok(ids.get('evidenceHead').textContent.includes('Ends'),'Legacy payload retains ends column');
  ids.get('motifRole').value='C';ids.get('motifRole').onchange();
  assert.equal(ids.get('motifExport').disabled,true,'No computable role gives empty motif');
  const v3=fixture();for(const t of v3.cleavage.termini)t.specificity='semi-specific';
  v3.analysis_files={pair_tests:'exports/pair_measurement_results.csv',pair_mapping_results:'exports/pair_results.csv',pair_exclusions:'exports/pair_excluded.csv',sample_ratios:'exports/neo_spanning_ratios.csv',pair_contrasts:'exports/neo_spanning_contrasts.csv',specificity:'exports/placement_specificity.csv',excluded_placements:'exports/excluded_placements.csv',source_input:'exports/source_input.xlsx'};
  v3.quantity_files={ratios:v3.analysis_files.sample_ratios};
  const v3ids=await run(v3);v3ids.get('sites').children[0].onclick();
  for(const file of Object.values(v3.analysis_files))assert.equal(v3ids.get('downloads').children.filter(a=>a.href===file).length,0,'Raw audit link kept outside primary UI: '+file);
  assert.ok(!v3ids.get('downloads').children.some(a=>a.href==='neo_spanning_contrasts.csv'),'No guessed analysis paths');
  assert.ok(v3ids.get('evidenceHead').textContent.includes('Specificity'));
  assert.ok(v3ids.get('evidenceHead').textContent.includes('Processing context'));
  const legacy=fixture();delete legacy.pair_view;delete legacy.merops;delete legacy.uniprot;
  const legacyIds=await run(legacy);assert.equal(legacyIds.get('motifExport').disabled,true);
  const grouped=fixture(),originalSite=grouped.cleavage.sites[0],duplicate={...originalSite,site_id:'alternate_mapping',record_id:99,accession:'ALT_MAPPING',bond_left:40,bond_right:41};
  grouped.cleavage.sites.push(duplicate);grouped.proteins.push({...grouped.proteins[0],record_id:99,accession:'ALT_MAPPING'});
  for(const key of ['evidence','termini'])grouped.cleavage[key].push(...grouped.cleavage[key].filter(e=>e.site_id===originalSite.site_id).map(e=>({...e,site_id:duplicate.site_id,record_id:99,accession:duplicate.accession,start:e.start+20,end:e.end+20})));
  grouped.matches.push(...grouped.matches.filter(e=>e.record_id===originalSite.record_id).map(e=>({...e,record_id:99,accession:duplicate.accession,start:e.start+20,end:e.end+20})));
  grouped.pair_view.ratios.push(...grouped.pair_view.ratios.filter(e=>e.site_id===originalSite.site_id).map(e=>({...e,site_id:duplicate.site_id,record_id:99,accession:duplicate.accession,bond_left:40,bond_right:41})));
  const groupedIds=await run(grouped),groupedRow=groupedIds.get('sites').children.find(r=>r.textContent.includes('ALT_MAPPING'));
  assert.ok(groupedRow.textContent.includes('FIXTURE_1'),'Both mappings listed in one row');
  assert.equal(groupedIds.get('sites').children.filter(r=>r.textContent.includes('ALT_MAPPING')||r.textContent.includes('FIXTURE_1')).length,1);
  groupedRow.onclick();assert.equal(groupedIds.get('mappingChoice').children.length,2);groupedIds.get('mappingChoice').value=duplicate.site_id;groupedIds.get('mappingChoice').onchange();assert.ok(groupedIds.get('bond').textContent.includes('ALT_MAPPING'));
  grouped.pair_view.statistics=[{site_id:originalSite.site_id,role:'N',neo_peptide_id:'mixed_a_neo',spanning_peptide_id:'mixed_a_span',comparison:'T_vs_C',family:'N_primary',status:'tested',p_value:.0123,adjusted_p_value:.0456}];
  const partialTests=await run(grouped);assert.ok(partialTests.get('sites').textContent.includes('Varies by mapping; see details'),'A test on one mapping must not be borrowed for another');
  grouped.pair_view.statistics.push({...grouped.pair_view.statistics[0],site_id:duplicate.site_id});
  const sameTests=await run(grouped);assert.ok(!sameTests.get('sites').textContent.includes('Varies by mapping; see details'));assert.ok(sameTests.get('sites').textContent.includes('0.0456'),'Identical attachment displayed once');
  grouped.pair_view.links=[...new Map(grouped.pair_view.ratios.map(r=>[JSON.stringify([r.site_id,r.role,r.neo_peptide_id,r.spanning_peptide_id]),r])).values()].map(r=>({...r,assignment_eligible:r.site_id!==duplicate.site_id,exclusion_reason:r.site_id===duplicate.site_id?'Outside shared source masters':''}));
  const gated=await run(grouped);assert.ok(!gated.get('sites').textContent.includes('ALT_MAPPING'),'Failed attribution outside primary overview');
  assert.ok(!gated.has('diagnosticSection'),'Mapping diagnostics is no longer rendered');
  assert.ok(grouped.pair_view.links.some(x=>!x.assignment_eligible),'Exclusion audit data is retained');
  console.log('Report in-memory DOM smoke: startup, motif scope, selectors, compact labels/audit, maps, specificity, explicit output links and legacy payload passed.');
  const args=process.argv.slice(2),currentTemplate=args.includes('--current-template');for(const file of args.filter(x=>x!=='--current-template'))await inspectFile(file,currentTemplate);
})().catch(error=>{console.error(error);process.exitCode=1});
