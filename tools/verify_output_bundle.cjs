// Read-only verification of a regenerated bundle; temporary Excel goes outside it.
'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict'),os=require('node:os'),{spawnSync}=require('node:child_process');
const Payload=require('../inst/templates/report_payload.js'),Workbook=require('../inst/templates/report_workbook.js');
const [source,newdir]=process.argv.slice(2);assert.ok(source&&newdir);
const hash=b=>crypto.createHash('sha256').update(b).digest('hex');
function read(file){const html=fs.readFileSync(file,'utf8');return {html,D:Payload.parseSync(html.match(/<script type="application\/json" id="data">([\s\S]*?)<\/script>/)[1])};}
(async()=>{
 const a=read(source),b=read(path.join(newdir,'report.html'));
 for(const key of ['peptides','matches','proteins','cleavage','merops','uniprot','site_summary','reference_access'])assert.deepEqual(b.D[key],a.D[key],key+' unchanged');
 for(const key of ['ratios','contrasts','measurements','samples','comparisons','links','statistics','excluded_statistics','inference_settings'])assert.deepEqual(b.D.pair_view?.[key],a.D.pair_view?.[key],key+' unchanged');
 assert.ok(!b.html.includes('About these results'));assert.ok(!b.html.includes('__ATELENTE_WORKBOOK_JS__'));assert.equal(b.D.excel_file,'atelente_results.xlsx');assert.equal(b.D.audit_file,'atelente_audit.zip');
 const sourceRoot=path.dirname(source);let csvs=0;
 const visit=(dir)=>{for(const item of fs.readdirSync(dir,{withFileTypes:true})){
   const file=path.join(dir,item.name);if(item.isDirectory()){visit(file);continue;}
   if(!item.isFile()||!item.name.endsWith('.csv'))continue;
   const relative=path.relative(sourceRoot,file),target=relative.startsWith('audit'+path.sep)?relative:path.join('audit',relative);
   assert.equal(hash(fs.readFileSync(file)),hash(fs.readFileSync(path.join(newdir,target))),relative+' source CSV unchanged');csvs++;
 }};visit(sourceRoot);
 const manifest=JSON.parse(fs.readFileSync(path.join(newdir,'run_manifest.json'),'utf8'));
 for(const [file,sha]of Object.entries(manifest.output_sha256))assert.equal(hash(fs.readFileSync(path.join(newdir,file))),sha);
 const index=JSON.parse(fs.readFileSync(path.join(newdir,'audit/audit_files.json'),'utf8'));for(const [file,sha]of Object.entries(index.sha256))assert.equal(hash(fs.readFileSync(path.join(newdir,'audit',file))),sha);
 const temp=fs.mkdtempSync(path.join(os.tmpdir(),'atelente-v6-verify-')),selected=path.join(temp,'selected.xlsx');
 fs.writeFileSync(selected,await Workbook.select(b.D.export_workbook,['Results','Peptides','Settings']));
 const script=path.join(temp,'verify.R'),main=path.resolve(newdir,'atelente_results.xlsx').replaceAll('\\','/');
 fs.writeFileSync(script,'a <- '+JSON.stringify(main)+'\nb <- '+JSON.stringify(selected.replaceAll('\\','/'))+'\nstopifnot(identical(readxl::excel_sheets(a),c("Results","Peptides","Settings")))\nfor(s in readxl::excel_sheets(a))stopifnot(identical(readxl::read_xlsx(a,sheet=s,trim_ws=FALSE),readxl::read_xlsx(b,sheet=s,trim_ws=FALSE)))\ncat("Real workbook matches offline selected download\\n")\n');
 const exe=process.env.RSCRIPT||(process.platform==='win32'?'C:/Program Files/R/R-4.3.0/bin/Rscript.exe':'Rscript'),r=spawnSync(exe,[script],{encoding:'utf8',maxBuffer:2e6});assert.equal(r.status,0,r.stdout+r.stderr);
 console.log(JSON.stringify({directory:path.resolve(newdir),unchanged_csvs:csvs,scientific_payload:'identical',audit_files:Object.keys(index.sha256).length,default_and_selected_excel:'identical readxl values/types',html_bytes:fs.statSync(path.join(newdir,'report.html')).size,excel_bytes:fs.statSync(main).size,zip_bytes:fs.statSync(path.join(newdir,'atelente_audit.zip')).size}));
})().catch(e=>{console.error(e.message);process.exitCode=1;});
