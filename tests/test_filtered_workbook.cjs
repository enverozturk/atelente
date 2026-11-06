'use strict';
// R creates ordinary typed XLSX inputs; readxl independently verifies browser output.
const assert=require('node:assert/strict'),fs=require('node:fs'),os=require('node:os'),path=require('node:path'),vm=require('node:vm');
const {spawnSync}=require('node:child_process'),W=require('../inst/templates/report_workbook.js');
const temp=fs.mkdtempSync(path.join(os.tmpdir(),'atelente-filtered-workbook-'));
const rscript=process.env.RSCRIPT||(process.platform==='win32'?'C:/Program Files/R/R-4.3.0/bin/Rscript.exe':'Rscript');
const root=path.resolve(__dirname,'..').replaceAll('\\','/'),dir=temp.replaceAll('\\','/');
function run(script,name){const file=path.join(temp,name+'.R');fs.writeFileSync(file,script);const r=spawnSync(rscript,[file],{encoding:'utf8'});assert.equal(r.status,0,r.stdout+'\n'+r.stderr);return r.stdout;}
const names=['Results','Peptides','Settings','All pairs','Other candidates','Annotations'];
const key=(cp,role,n,s)=>JSON.stringify([cp,role,n,s]);
class Element{
  constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.textContent='';}
  append(...x){this.children.push(...x);}replaceChildren(...x){this.children=x;}setAttribute(k,v){this.attributes[k]=v;}
  click(){this.onclick?.();}
}
(async()=>{
  run(`source(${JSON.stringify(root+'/R/excel_locale.R')});source(${JSON.stringify(root+'/R/excel.R')})
  overview <- data.frame(Comparison=c('KO_vs_WT','KO_vs_WT','B_vs_A'),Protein=c('P1; P2','P3','P4'),Role='Neo-N',
    'Cleavage window'=c('AAAA | BBBB','CCCC | DDDD','EEEE | FFFF'),
    'Neo log2FC'=c(1.23456789012345,3.331e-16,-1), 'Spanning log2FC'=c(0.2,0,0),
    'Ratio log2FC'=c(1.03456789012345,3.331e-16,-1),Pattern='DEFAULT_PATTERN',
    'Control measured'=c(3L,1L,2L),'Treatment measured'=c(5L,1L,2L),
    'P-value'=c(1.23456789012345e-120,NA,0.05),'BH-adjusted P'=c(2.34567890123456e-120,NA,0.1),
    'Neo form ID'=c('00123','N2','N3'),'Spanning form ID'=c('S1','S2','S3'),
    Control=c('WT','WT','A'),Treatment=c('KO','KO','B'),'Control total'=c(6L,6L,2L),'Treatment total'=c(5L,5L,2L),
    'Site IDs'=c('site1; site2','site3','site4'),'Export key'=c(${JSON.stringify(key('KO_vs_WT','N','00123','S1'))},${JSON.stringify(key('KO_vs_WT','N','N2','S2'))},${JSON.stringify(key('B_vs_A','N','N3','S3'))}),
    'Passes default coverage'=c(TRUE,FALSE,TRUE),check.names=FALSE)
  pairs <- overview[c(1,1,1,1,1,2,3),,drop=FALSE]
  pairs[['Site IDs']] <- c('site1','site2','site1','site1','site1','site3','site4')
  pairs[['Neo form ID']][3:5] <- c('ALT','WRONG_ROLE','WRONG_COMPARISON')
  pairs[['Role']][4] <- 'Neo-C';pairs[['Comparison']][5] <- 'B_vs_A'
  pairs[['Ratio log2FC']][3] <- -1.125
  ids <- c('00123','S1','ALT','N2','S2','N3','S3','WRONG_ROLE','WRONG_COMPARISON','UNRELATED')
  peps <- data.frame('Form ID'=ids,Modifications=c('alpha & beta','',rep(NA,8)), 'Quantity: sample 1'=c(0,NA,1.23456789012345,1,2,3,4,5,6,7),check.names=FALSE)
  settings <- data.frame(Section=c('Analysis','Comparison: KO_vs_WT','Export','Export'),Setting=c('Frozen selection k','Default display coverage k','Scope','Ratio threshold'),Value=c('3','3','All results','1'))
  tables <- list(Results=overview[c(1,3),1:14],Peptides=peps,Settings=settings,'All pairs'=pairs,'Other candidates'=peps[10,,drop=FALSE],
    Annotations=data.frame(site_id=c('site1','site2','site3','site4'),Description=c('ref1','ref2','ref3','ref4')), '__Results'=overview)
  pm_write_excel_tables(tables,${JSON.stringify(dir+'/source.xlsx')})`, 'create');
  const workbook={sheets:names,internal_sheets:['__Results'],base64:fs.readFileSync(path.join(temp,'source.xlsx')).toString('base64')};
  const view={rows:[{key:key('KO_vs_WT','N','N2','S2'),pattern:'Below ratio threshold'},{key:key('KO_vs_WT','N','00123','S1'),pattern:'Ratio increased'}],
    settings:{comparison:'KO_vs_WT',control:'WT',treatment:'KO',controlTotal:6,treatmentTotal:5,threshold:0.5,coverage:1,analysisSelectionK:3,search:'P1 | P3',pattern:'All patterns',rowCount:2},
    pairPatterns:[{key:key('KO_vs_WT','N','00123','S1'),pattern:'Ratio increased'},{key:key('KO_vs_WT','N','ALT','S1'),pattern:'Ratio decreased'},{key:key('KO_vs_WT','N','N2','S2'),pattern:'Below ratio threshold'}]};
  fs.writeFileSync(path.join(temp,'filtered.xlsx'),await W.select(workbook,names,{view}));
  fs.writeFileSync(path.join(temp,'mixed.xlsx'),await W.select(workbook,['Results'],{view:{...view,rows:[view.rows[1],{key:key('B_vs_A','N','N3','S3'),pattern:'Ratio decreased'}]}}));
  fs.writeFileSync(path.join(temp,'empty.xlsx'),await W.select(workbook,['Results'],{view:{...view,rows:[]}}));
  const chosen=await W.select(workbook,['Results']);assert.ok(!Buffer.from(chosen).toString('utf8').includes('__Results'));
  await assert.rejects(W.select(workbook,['__Results']),/Select at least/);
  await assert.rejects(W.select(workbook,['Results'],{view:{...view,rows:[view.rows[0],view.rows[0]]}}),/duplicate/);
  await assert.rejects(W.select(workbook,['Results'],{view:{...view,rows:[{key:'bad',pattern:'bad'}]}}),/does not match/);
  run(`a <- ${JSON.stringify(dir+'/source.xlsx')};b <- ${JSON.stringify(dir+'/filtered.xlsx')}
  read <- function(f,s) as.data.frame(readxl::read_xlsx(f,sheet=s,trim_ws=FALSE))
  x <- read(b,'Results');o <- read(a,'__Results')
  stopifnot(nrow(x)==2L,identical(x[['Neo form ID']],c('N2','00123')),
    is.numeric(x[['WT measured (of 6)']]),identical(x[['WT measured (of 6)']],c(1,3)),
    identical(x[['KO measured (of 5)']],c(1,5)),
    identical(x[['Ratio log2FC']],o[['Ratio log2FC']][c(2,1)]),
    identical(x[['P-value']],o[['P-value']][c(2,1)]),
    identical(x[['BH-adjusted P']],o[['BH-adjusted P']][c(2,1)]),
    identical(x[['Pattern']],c('Below ratio threshold','Ratio increased')),
    !any(c('Export key','Site IDs','Control','Treatment','Control total','Treatment total','Passes default coverage') %in% names(x)))
  p <- read(b,'Peptides');stopifnot(setequal(p[['Form ID']],c('00123','S1','ALT','N2','S2')),
    p[['Quantity: sample 1']][p[['Form ID']]=='00123']==0,is.na(p[['Quantity: sample 1']][p[['Form ID']]=='S1']))
  q <- read(b,'All pairs');stopifnot(nrow(q)==4L,all(q$Role=='Neo-N'),all(q$Comparison=='KO_vs_WT'),
    !any(q[['Neo form ID']] %in% c('WRONG_ROLE','WRONG_COMPARISON')),
    q$Pattern[q[['Neo form ID']]=='ALT']=='Ratio decreased')
  z <- read(b,'Annotations');stopifnot(setequal(z$site_id,c('site1','site2','site3')))
  t <- read(b,'Settings');stopifnot(t$Value[t$Section=='Export' & t$Setting=='Ratio threshold']=='0.5',
    t$Value[t$Section=='Export' & t$Setting=='Display coverage k']=='1',
    t$Value[t$Section=='Analysis' & t$Setting=='Frozen selection k']=='3',
    !any(t$Value[t$Section=='Export']=='All results'),t$Value[t$Section=='Export' & t$Setting=='Search']=='P1 | P3')
  m <- ${JSON.stringify(dir+'/mixed.xlsx')};stopifnot(identical(readxl::excel_sheets(m),c('Results','Settings')))
  y <- read(m,'Results');stopifnot(all(c('Comparison','Control measured','Treatment measured') %in% names(y)),nrow(y)==2L)
  e <- ${JSON.stringify(dir+'/empty.xlsx')};stopifnot(nrow(read(e,'Results'))==0L,identical(readxl::excel_sheets(e),c('Results','Settings')))
  cat('Filtered Excel: saved numeric precision, below-default rows, ordering, typed coverage, metadata, alternate scope and missing/zero passed\\n')`, 'verify');
  // UI state is captured at click time and the standalone browser UMD works.
  const ids=new Map(['downloads','tableChooser','tableChoices','downloadChosen','downloadStatus'].map(id=>[id,new Element('div')]));
  const doc={getElementById:id=>ids.get(id),createElement:tag=>new Element(tag)};
  let downloaded,viewCalls=0;
  const context={TextEncoder,TextDecoder,Uint8Array,Uint32Array,DataView,DecompressionStream,Blob,Response,atob,
    URL:{createObjectURL(blob){downloaded=blob;return 'blob:local';},revokeObjectURL(){}},setTimeout:fn=>fn()};
  vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../inst/templates/report_workbook.js'),'utf8'),context);
  const controls=context.ReportWorkbook.mount({excel_file:'atelente_results.xlsx',audit_file:'atelente_audit.zip',export_workbook:workbook},doc,{getView:()=>{viewCalls++;return view;}});
  assert.equal(ids.get('downloads').children.length,3);
  await controls.downloadCurrent();assert.equal(viewCalls,1);assert.match(ids.get('downloadStatus').textContent,/2 results/);
  assert.deepEqual(Buffer.from(await downloaded.arrayBuffer()),Buffer.from(await W.select(workbook,['Results','Peptides','Settings'],{view})));
  const scope=ids.get('tableChoices').children[0].children[0];assert.equal(scope.id,'exportScope');scope.value='current';
  const checks=ids.get('tableChoices').children.slice(1).map(label=>label.children[0]);checks.forEach(c=>c.checked=c.value==='Results');
  await ids.get('downloadChosen').onclick();assert.equal(viewCalls,2);assert.match(ids.get('downloadStatus').textContent,/2 tables/);
  console.log('Filtered workbook UI scope, snapshot capture and current-view download passed.');
})().catch(error=>{console.error(error);process.exitCode=1;});
