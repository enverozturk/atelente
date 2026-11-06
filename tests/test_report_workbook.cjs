'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),os=require('node:os'),path=require('node:path'),zlib=require('node:zlib'),vm=require('node:vm');
const {spawnSync}=require('node:child_process'),W=require('../inst/templates/report_workbook.js');
const temp=fs.mkdtempSync(path.join(os.tmpdir(),'atelente-workbook-'));
const REL='http://schemas.openxmlformats.org/officeDocument/2006/relationships',MAIN='http://schemas.openxmlformats.org/spreadsheetml/2006/main';
const crcTable=Uint32Array.from({length:256},(_,c)=>{for(let i=0;i<8;i++)c=c&1?0xedb88320^(c>>>1):c>>>1;return c>>>0;});
function crc(bytes){let c=0xffffffff;for(const b of bytes)c=crcTable[(c^b)&255]^(c>>>8);return (c^0xffffffff)>>>0;}
function zip(parts,method=8){
  const local=[],central=[];let offset=0;
  for(const [name,value] of parts){const n=Buffer.from(name),plain=Buffer.from(value),data=method===8?zlib.deflateRawSync(plain):plain,c=crc(plain),h=Buffer.alloc(30),d=Buffer.alloc(46);
    h.writeUInt32LE(0x04034b50);h.writeUInt16LE(20,4);h.writeUInt16LE(method,8);h.writeUInt32LE(c,14);h.writeUInt32LE(data.length,18);h.writeUInt32LE(plain.length,22);h.writeUInt16LE(n.length,26);
    d.writeUInt32LE(0x02014b50);d.writeUInt16LE(20,4);d.writeUInt16LE(20,6);d.writeUInt16LE(method,10);d.writeUInt32LE(c,16);d.writeUInt32LE(data.length,20);d.writeUInt32LE(plain.length,24);d.writeUInt16LE(n.length,28);d.writeUInt32LE(offset,42);
    local.push(h,n,data);central.push(d,n);offset+=h.length+n.length+data.length;
  }
  const cd=Buffer.concat(central),end=Buffer.alloc(22);end.writeUInt32LE(0x06054b50);end.writeUInt16LE(parts.length,8);end.writeUInt16LE(parts.length,10);end.writeUInt32LE(cd.length,12);end.writeUInt32LE(offset,16);
  return Buffer.concat([...local,cd,end]);
}
function unzip(bytes){const b=Buffer.from(bytes),end=b.length-22,files=new Map();assert.equal(b.readUInt32LE(end),0x06054b50);
  let cursor=b.readUInt32LE(end+16);for(let i=0;i<b.readUInt16LE(end+10);i++){
    assert.equal(b.readUInt32LE(cursor),0x02014b50);const n=b.readUInt16LE(cursor+28),extra=b.readUInt16LE(cursor+30),comment=b.readUInt16LE(cursor+32),name=b.toString('utf8',cursor+46,cursor+46+n),loc=b.readUInt32LE(cursor+42),size=b.readUInt32LE(cursor+20),method=b.readUInt16LE(cursor+10);
    const start=loc+30+b.readUInt16LE(loc+26)+b.readUInt16LE(loc+28),compressed=b.subarray(start,start+size),data=method===8?zlib.inflateRawSync(compressed):compressed;
    assert.equal(data.length,b.readUInt32LE(cursor+24));assert.equal(crc(data),b.readUInt32LE(cursor+16));assert.ok(!files.has(name));files.set(name,data.toString('utf8'));cursor+=46+n+extra+comment;
  }return files;
}
const names=['Results','Peptides','Settings','Other candidates','All pairs','Annotations'];
const shared=['ID','Value','Flag','00123','\u03b1 &amp; \u03b2 &lt;literal&gt;','=SUM(1,2)','DO_NOT_EXPORT_UNSELECTED',' secret leading space '];
const sheet=(secret=false)=>'<worksheet xmlns="'+MAIN+'" xmlns:r="'+REL+'"><sheetViews><sheetView workbookViewId="0" tabSelected="true"><pane ySplit="1" xSplit="1" topLeftCell="B2" state="frozen"/></sheetView></sheetViews><cols><col min="1" max="3" width="24" customWidth="1"/></cols><sheetData><row r="1"><c r="A1" s="1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c><c r="C1" t="s"><v>2</v></c></row><row r="2"><c r="A2" s="2" t="s"><v>'+(secret?6:3)+'</v></c><c r="B2" t="n"><v>1.23456789012345e-120</v></c><c r="C2" t="b"><v>1</v></c></row><row r="3"><c r="A3" t="s"><v>4</v></c><c r="B3"/><c r="C3" t="b"><v>0</v></c></row><row r="4"><c r="A4" t="s"><v>5</v></c><c r="B4" t="n"><v>0</v></c><c r="C4" t="s"><v>7</v></c></row></sheetData><autoFilter ref="A1:C4"/><pageSetup orientation="portrait" r:id="rId1"/></worksheet>';
const parts=[['xl/workbook.xml','<workbook xmlns="'+MAIN+'" xmlns:r="'+REL+'"><workbookPr date1904="true"/><sheets>'+names.map((n,i)=>'<sheet name="'+n+'" sheetId="'+(i+1)+'" r:id="rId'+(i+1)+'"/>').join('')+'</sheets></workbook>'],
 ['xl/_rels/workbook.xml.rels','<Relationships>'+names.map((_,i)=>'<Relationship Id="rId'+(i+1)+'" Type="'+REL+'/worksheet" Target="worksheets/sheet'+(i+1)+'.xml"/>').join('')+'<Relationship Id="strings" Type="'+REL+'/sharedStrings" Target="sharedStrings.xml"/><Relationship Id="styles" Type="'+REL+'/styles" Target="styles.xml"/><Relationship Id="theme" Type="'+REL+'/theme" Target="theme/theme1.xml"/></Relationships>'],
 ['xl/sharedStrings.xml','<sst xmlns="'+MAIN+'">'+shared.map(s=>'<si><t xml:space="preserve">'+s+'</t></si>').join('')+'</sst>'],
 ['xl/styles.xml','<styleSheet xmlns="'+MAIN+'"><cellXfs count="3"><xf/><xf/><xf numFmtId="49"/></cellXfs></styleSheet>'],
 ['xl/theme/theme1.xml','<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="Atelente"/>'],
 ...names.map((_,i)=>['xl/worksheets/sheet'+(i+1)+'.xml',sheet(i===1)])];
const fixture={sheets:names,base64:zip(parts).toString('base64')};
class Element{
  constructor(tag){this.tagName=tag;this.children=[];this.attributes={};this.textContent='';}
  append(...x){this.children.push(...x);}replaceChildren(...x){this.children=x;}setAttribute(k,v){this.attributes[k]=v;}
  click(){this.onclick?.();}
}
(async()=>{
  const chosen=await W.select(fixture,['Results','Settings']),files=unzip(chosen);
  assert.equal(files.size,8);assert.ok(!files.has('xl/sharedStrings.xml'));
  assert.ok(!Buffer.from(chosen).toString('utf8').includes('DO_NOT_EXPORT_UNSELECTED'),'Unselected shared strings do not leak');
  assert.ok(!files.get('xl/workbook.xml').includes('Peptides'));assert.ok(files.get('xl/workbook.xml').includes('date1904="true"'));
  const xml=files.get('xl/worksheets/sheet1.xml');
  for(const text of ['t="inlineStr"','xml:space="preserve">00123','1.23456789012345e-120','t="b"><v>1','t="b"><v>0','<c r="B3"/>','<v>0</v>','\u03b1 &amp; \u03b2 &lt;literal&gt;','=SUM(1,2)',' secret leading space ','state="frozen"','<autoFilter ref="A1:C4"'])assert.ok(xml.includes(text),text);
  assert.ok(!xml.includes('r:id='));assert.ok(!xml.includes('<f>'));assert.ok(files.get('xl/workbook.xml').includes('$A$1:$C$4'));
  assert.ok(files.get('xl/worksheets/sheet2.xml').includes('tabSelected="false"'));assert.equal(files.get('xl/styles.xml'),parts.find(x=>x[0]==='xl/styles.xml')[1]);
  assert.deepEqual(Buffer.from(await W.select(fixture,['Results','Settings'])),Buffer.from(chosen),'Deterministic offline output');
  const every=unzip(await W.select(fixture,names));assert.ok([...every.values()].some(x=>x.includes('DO_NOT_EXPORT_UNSELECTED')));
  for(const selection of [[],['Unknown'],['Results','Results']])await assert.rejects(W.select(fixture,selection),/Select at least/);
  await assert.rejects(W.select({...fixture,sheets:['Results']},['Results']),/table list/);
  await assert.rejects(W.select({...fixture,base64:'AAAA'},['Results']),/ZIP/);
  await assert.rejects(W.select({...fixture,base64:zip([...parts,['../escape.xml','bad']]).toString('base64')},['Results']),/Unsafe/);
  await assert.rejects(W.select({...fixture,base64:zip([...parts,['xl/workbook.xml','duplicate']]).toString('base64')},['Results']),/Duplicate/);
  const corrupt=zip(parts,0),offset=corrupt.indexOf(Buffer.from('date1904="true"'));corrupt[offset]='X'.charCodeAt(0);
  await assert.rejects(W.select({...fixture,base64:corrupt.toString('base64')},['Results']),/checksum/);
  const entity=parts.map(([p,v])=>[p,p==='xl/workbook.xml'?'<!DOCTYPE root [<!ENTITY secret SYSTEM "file:///secret">]>'+v:v]);
  await assert.rejects(W.select({...fixture,base64:zip(entity).toString('base64')},['Results']),/declarations/);
  const badIndex=parts.map(([p,v])=>[p,p==='xl/worksheets/sheet1.xml'?v.replace('<v>3</v>','<v>99999</v>'):v]);
  await assert.rejects(W.select({...fixture,base64:zip(badIndex).toString('base64')},['Results']),/shared-string/);
  const afterBlank=parts.map(([p,v])=>[p,p==='xl/worksheets/sheet1.xml'?v.replace('<c r="C3" t="b"><v>0</v></c>','<c r="C3" t="s"><v>4</v></c>'):v]);
  const converted=unzip(await W.select({...fixture,base64:zip(afterBlank).toString('base64')},['Results'])).get('xl/worksheets/sheet1.xml');
  assert.ok(converted.includes('<c r="B3"/><c r="C3" t="inlineStr"><is>'),'Empty cells do not consume the next shared-string cell');
  const emptyString=parts.map(([p,v])=>[p,p==='xl/sharedStrings.xml'?v.replace('<si><t xml:space="preserve">ID</t></si>','<si/>'):v]);
  const emptyConverted=unzip(await W.select({...fixture,base64:zip(emptyString).toString('base64')},['Results'])).get('xl/worksheets/sheet1.xml');
  assert.ok(emptyConverted.includes('<c r="A1" s="1" t="inlineStr"><is></is>'),'Empty shared-string entries retain their index');
  assert.ok(emptyConverted.includes('00123'),'Later string indices remain stable');
  const formula=parts.map(([p,v])=>[p,p==='xl/worksheets/sheet1.xml'?v.replace('<v>3</v>','<f>Peptides!A1</f><v>3</v>'):v]);
  await assert.rejects(W.select({...fixture,base64:zip(formula).toString('base64')},['Results']),/ordinary/);

  const ids=new Map(['downloads','tableChooser','tableChoices','downloadChosen','downloadStatus'].map(id=>[id,new Element('div')]));
  const doc={getElementById:id=>ids.get(id),createElement:tag=>new Element(tag)};
  W.mount({excel_file:'atelente_results.xlsx',audit_file:'atelente_audit.zip',export_workbook:fixture},doc);
  assert.equal(ids.get('downloads').children.length,3);assert.equal(ids.get('downloads').children[0].href,'atelente_results.xlsx');assert.equal(ids.get('downloads').children[2].href,'atelente_audit.zip');
  assert.equal(ids.get('tableChooser').hidden,true);ids.get('downloads').children[1].click();assert.equal(ids.get('tableChooser').hidden,false);
  const checks=ids.get('tableChoices').children.map(x=>x.children[0]);assert.deepEqual(checks.filter(c=>c.checked).map(c=>c.value),['Results','Peptides','Settings']);
  checks.forEach(c=>c.checked=false);checks[0].onchange();assert.equal(ids.get('downloadChosen').disabled,true);
  checks[0].checked=true;checks[0].onchange();assert.equal(ids.get('downloadChosen').disabled,false);
  W.mount({},doc);assert.equal(ids.get('downloads').children[0].disabled,true);
  assert.throws(()=>W.mount({excel_file:'https://example.com/private.xlsx'},doc),/Unsafe/);
  // Browser UMD execution has no Node require dependency.
  let downloaded,revoked=false;
  const context={TextEncoder,TextDecoder,Uint8Array,Uint32Array,DataView,DecompressionStream,Blob,Response,atob,
    URL:{createObjectURL(blob){downloaded=blob;return 'blob:local-export';},revokeObjectURL(url){assert.equal(url,'blob:local-export');revoked=true;}},setTimeout:fn=>fn()};
  vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../inst/templates/report_workbook.js'),'utf8'),context);
  assert.deepEqual(Buffer.from(await context.ReportWorkbook.select(fixture,['Results'])),Buffer.from(await W.select(fixture,['Results'])));
  context.ReportWorkbook.mount({excel_file:'atelente_results.xlsx',audit_file:'atelente_audit.zip',export_workbook:fixture},doc);
  ids.get('tableChoices').children.forEach(label=>label.children[0].checked=label.children[0].value==='Results');
  await ids.get('downloadChosen').onclick();
  assert.equal(ids.get('downloadStatus').textContent,'Excel prepared · 1 tables');assert.equal(ids.get('downloadChosen').disabled,false);assert.equal(revoked,true);
  assert.deepEqual(Buffer.from(await downloaded.arrayBuffer()),Buffer.from(await W.select(fixture,['Results'])));

  // Produce a real openxlsx workbook, then verify selection with the independent
  // readxl reader. openxlsx 4.2.5's inline-string reader drops XML attributes poorly;
  // preserving valid OOXML (including xml:space) is required for Excel/readxl.
  const rscript=process.env.RSCRIPT||(process.platform==='win32'?'C:/Program Files/R/R-4.3.0/bin/Rscript.exe':'Rscript');
  const probe=spawnSync(rscript,['--version'],{encoding:'utf8'});
  if(probe.error?.code==='ENOENT')console.log('R roundtrip skipped: Rscript not installed.');
  else{
    const root=path.resolve(__dirname,'..').replaceAll('\\','/'),dir=temp.replaceAll('\\','/'),create=path.join(temp,'create.R'),verify=path.join(temp,'verify.R');
    fs.writeFileSync(create,'source('+JSON.stringify(root+'/R/excel_locale.R')+')\nsource('+JSON.stringify(root+'/R/excel.R')+')\ntables <- list(Results=data.frame(ID=c("00123","\\u03b1 & \\u03b2 <literal>","=SUM(1,2)"),Value=c(1.23456789012345e-120,NA,0),Flag=c(TRUE,FALSE,NA),check.names=FALSE),Peptides=data.frame(Secret="DO_NOT_EXPORT_UNSELECTED"),Settings=data.frame(Key="space",Value=" padded "))\npm_write_excel_tables(tables,'+JSON.stringify(dir+'/source.xlsx')+')\n');
    const generated=spawnSync(rscript,[create],{encoding:'utf8'});assert.equal(generated.status,0,generated.stderr);
    const original=fs.readFileSync(path.join(temp,'source.xlsx')),selected=await W.select({sheets:['Results','Peptides','Settings'],base64:original.toString('base64')},['Results','Settings']);
    fs.writeFileSync(path.join(temp,'selected.xlsx'),selected);assert.ok(!Buffer.from(selected).toString('utf8').includes('DO_NOT_EXPORT_UNSELECTED'));
    fs.writeFileSync(verify,'a <- '+JSON.stringify(dir+'/source.xlsx')+'\nb <- '+JSON.stringify(dir+'/selected.xlsx')+'\nstopifnot(identical(readxl::excel_sheets(b),c("Results","Settings")))\nfor(s in c("Results","Settings")) stopifnot(identical(readxl::read_xlsx(a,sheet=s,trim_ws=FALSE),readxl::read_xlsx(b,sheet=s,trim_ws=FALSE)))\ncat("openxlsx source -> offline selection -> readxl roundtrip passed\\n")\n');
    const checked=spawnSync(rscript,[verify],{encoding:'utf8'});assert.equal(checked.status,0,checked.stdout+'\n'+checked.stderr);process.stdout.write(checked.stdout);
  }
  console.log('Workbook selection: offline ZIP/CRC, no unselected-string leakage, OOXML data/styles/freeze/filter preservation, table controls and malformed-input rejection passed.');
})().catch(error=>{console.error(error);process.exitCode=1;});
