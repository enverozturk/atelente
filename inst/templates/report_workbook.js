/* Offline selection of ordinary Atelente/openxlsx worksheets. No network or formulas. */
(function(root,factory){
  const api=factory();
  if(typeof module==='object'&&module.exports)module.exports=api;else root.ReportWorkbook=api;
})(typeof globalThis!=='undefined'?globalThis:this,function(){
  'use strict';
  const encoder=new TextEncoder(),decoder=new TextDecoder('utf-8',{fatal:true});
  const LIMIT=512*1024*1024,MAIN='http://schemas.openxmlformats.org/spreadsheetml/2006/main',
    REL='http://schemas.openxmlformats.org/officeDocument/2006/relationships',PKG='http://schemas.openxmlformats.org/package/2006/relationships';
  const crcTable=Uint32Array.from({length:256},(_,n)=>{for(let k=0;k<8;k++)n=n&1?0xedb88320^(n>>>1):n>>>1;return n>>>0;});
  function crc32(bytes){let crc=0xffffffff;for(const byte of bytes)crc=crcTable[(crc^byte)&255]^(crc>>>8);return (crc^0xffffffff)>>>0;}
  function safePath(name){
    if(typeof name!=='string'||!name||/^[\\/]|[\\:\u0000-\u001f]|%/u.test(name)||name.split('/').some(p=>p==='..'||p==='.'||!p))throw Error('Unsafe workbook ZIP path.');
    return name;
  }
  function base64(value){
    if(typeof value==='string')value=value.replace(/[\r\n\t ]/g,'');
    if(typeof value!=='string'||!value.length||value.length%4||!/^[A-Za-z0-9+/]*={0,2}$/.test(value)||value.length>LIMIT*4/3)throw Error('Invalid workbook payload.');
    return Uint8Array.from(atob(value),c=>c.charCodeAt(0));
  }
  function archive(bytes){
    const view=new DataView(bytes.buffer,bytes.byteOffset,bytes.byteLength),size=bytes.length;
    const inside=(offset,length)=>Number.isSafeInteger(offset)&&offset>=0&&length>=0&&offset+length<=size;
    const u16=o=>view.getUint16(o,true),u32=o=>view.getUint32(o,true);
    let end=-1;for(let p=size-22;p>=Math.max(0,size-65557);p--)if(u32(p)===0x06054b50&&p+22+u16(p+20)===size){end=p;break;}
    if(end<0||u16(end+4)||u16(end+6)||u16(end+8)!==u16(end+10))throw Error('Unsupported or incomplete workbook ZIP.');
    const count=u16(end+10),length=u32(end+12),start=u32(end+16);
    if(count===65535||start===0xffffffff||length===0xffffffff||!inside(start,length)||start+length!==end)throw Error('Unsupported workbook ZIP directory.');
    const entries=new Map();let cursor=start,total=0;
    for(let i=0;i<count;i++){
      if(!inside(cursor,46)||u32(cursor)!==0x02014b50)throw Error('Invalid workbook ZIP directory entry.');
      const flags=u16(cursor+8),method=u16(cursor+10),crc=u32(cursor+16),compressed=u32(cursor+20),uncompressed=u32(cursor+24),
        n=u16(cursor+28),extra=u16(cursor+30),comment=u16(cursor+32),disk=u16(cursor+34),local=u32(cursor+42);
      if(!inside(cursor+46,n+extra+comment)||cursor+46+n+extra+comment>start+length||disk||flags&1||flags&64||![0,8].includes(method)||
        compressed===0xffffffff||uncompressed===0xffffffff||(total+=uncompressed)>LIMIT)throw Error('Unsupported workbook ZIP entry.');
      const name=safePath(decoder.decode(bytes.subarray(cursor+46,cursor+46+n)));
      if(entries.has(name))throw Error('Duplicate workbook ZIP path.');
      if(!inside(local,30)||local>=start||u32(local)!==0x04034b50||u16(local+6)!==flags||u16(local+8)!==method)throw Error('Invalid workbook ZIP local entry.');
      const localName=u16(local+26),localExtra=u16(local+28),data=local+30+localName+localExtra;
      if(!inside(local+30,localName+localExtra)||!inside(data,compressed)||data+compressed>start||
        decoder.decode(bytes.subarray(local+30,local+30+localName))!==name)throw Error('Workbook ZIP entry boundary mismatch.');
      if(!(flags&8)&&(u32(local+14)!==crc||u32(local+18)!==compressed||u32(local+22)!==uncompressed))throw Error('Workbook ZIP entry metadata mismatch.');
      entries.set(name,{method,crc,compressed,uncompressed,data});cursor+=46+n+extra+comment;
    }
    if(cursor!==start+length)throw Error('Workbook ZIP directory length mismatch.');
    const cache=new Map();
    async function read(name){
      if(cache.has(name))return cache.get(name);
      const entry=entries.get(name);if(!entry)throw Error('Missing workbook part: '+name);
      const source=bytes.subarray(entry.data,entry.data+entry.compressed);let output;
      if(entry.method===0)output=source;
      else{
        if(typeof DecompressionStream!=='function')throw Error('Choose tables requires a current browser with offline decompression support.');
        const reader=new Blob([source]).stream().pipeThrough(new DecompressionStream('deflate-raw')).getReader(),chunks=[];let length=0;
        try{for(;;){const {value,done}=await reader.read();if(done)break;length+=value.length;
          if(length>entry.uncompressed||length>LIMIT){await reader.cancel();throw Error('Workbook decompression limit exceeded.');}chunks.push(value);}}
        finally{reader.releaseLock();}
        output=new Uint8Array(length);let offset=0;for(const chunk of chunks){output.set(chunk,offset);offset+=chunk.length;}
      }
      if(output.length!==entry.uncompressed||crc32(output)!==entry.crc)throw Error('Workbook ZIP checksum mismatch.');
      cache.set(name,output);return output;
    }
    return {entries,read};
  }
  function xmlText(bytes){const text=decoder.decode(bytes);if(/<!DOCTYPE|<!ENTITY/i.test(text))throw Error('Unsupported workbook XML declarations.');return text;}
  function decodeXml(text){return text.replace(/&(#x[0-9a-f]+|#[0-9]+|amp|lt|gt|quot|apos);/gi,(_,x)=>
    x[0]==='#'?String.fromCodePoint(x[1].toLowerCase()==='x'?parseInt(x.slice(2),16):Number(x.slice(1))):({amp:'&',lt:'<',gt:'>',quot:'"',apos:"'"})[x]);}
  const escapeXml=text=>String(text).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'})[c]);
  function attrs(text){const out={};for(const match of text.matchAll(/([A-Za-z_][\w.:-]*)\s*=\s*(["'])(.*?)\2/gs))out[match[1]]=decodeXml(match[3]);return out;}
  function resolvePart(base,target){
    if(!target||/[\\:\u0000-\u001f]|%/.test(target))throw Error('Unsafe workbook relationship.');
    const parts=target[0]==='/'?[]:base.split('/').slice(0,-1);
    for(const part of target.replace(/^\//,'').split('/')){if(part==='..'){if(!parts.length)throw Error('Unsafe workbook relationship.');parts.pop();}else if(part&&part!=='.')parts.push(part);}
    return safePath(parts.join('/'));
  }
  function relationships(text,base){const map=new Map();for(const m of text.matchAll(/<Relationship\b([^>]*?)\/?\s*>/g)){
    const a=attrs(m[1]);if(!a.Id||map.has(a.Id)||!a.Type||a.TargetMode==='External')throw Error('Unsupported workbook relationship.');
    map.set(a.Id,{type:a.Type,path:resolvePart(base,a.Target)});
  }return map;}
  function storeZip(parts){
    let localSize=0,centralSize=0;const encoded=parts.map(([name,value])=>{safePath(name);const filename=encoder.encode(name),bytes=typeof value==='string'?encoder.encode(value):value;
      localSize+=30+filename.length+bytes.length;centralSize+=46+filename.length;return {filename,bytes,crc:crc32(bytes)};});
    if(encoded.length>=65535||localSize+centralSize+22>LIMIT)throw Error('Selected workbook exceeds offline export size limit.');
    const out=new Uint8Array(localSize+centralSize+22),view=new DataView(out.buffer);let cursor=0;
    const u16=(o,v)=>view.setUint16(o,v,true),u32=(o,v)=>view.setUint32(o,v,true);
    for(const e of encoded){e.offset=cursor;u32(cursor,0x04034b50);u16(cursor+4,20);u16(cursor+6,0x800);u16(cursor+12,33);
      u32(cursor+14,e.crc);u32(cursor+18,e.bytes.length);u32(cursor+22,e.bytes.length);u16(cursor+26,e.filename.length);
      out.set(e.filename,cursor+30);out.set(e.bytes,cursor+30+e.filename.length);cursor+=30+e.filename.length+e.bytes.length;}
    for(const e of encoded){u32(cursor,0x02014b50);u16(cursor+4,20);u16(cursor+6,20);u16(cursor+8,0x800);u16(cursor+14,33);
      u32(cursor+16,e.crc);u32(cursor+20,e.bytes.length);u32(cursor+24,e.bytes.length);u16(cursor+28,e.filename.length);u32(cursor+42,e.offset);
      out.set(e.filename,cursor+46);cursor+=46+e.filename.length;}
    u32(cursor,0x06054b50);u16(cursor+8,encoded.length);u16(cursor+10,encoded.length);u32(cursor+12,centralSize);u32(cursor+16,localSize);
    return out;
  }
  // Parse ordinary data cells without converting numeric XML payloads to JS
  // numbers. Filtering/reordering therefore preserves every saved numeric digit.
  const columnIndex=letters=>{let n=0;for(const c of letters)n=n*26+c.charCodeAt(0)-64;return n-1;};
  const columnName=index=>{let s='';for(let n=index+1;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
  function cellValue(cell){
    if(!cell)return '';
    if(cell.attributes.t==='inlineStr')return [...cell.body.matchAll(/<t\b[^>]*>([\s\S]*?)<\/t>/g)].map(m=>decodeXml(m[1])).join('');
    const m=cell.body.match(/<v\b[^>]*>([\s\S]*?)<\/v>/);return m?decodeXml(m[1]):'';
  }
  function textCell(value,original){return {attributes:{...(original?.attributes||{}),t:'inlineStr'},body:'<is><t xml:space="preserve">'+escapeXml(value??'')+'</t></is>'};}
  function dataTable(sheet){
    const data=sheet.match(/<sheetData\b[^>]*>([\s\S]*?)<\/sheetData>/);
    if(!data)throw Error('Missing workbook table data.');
    const rows=[...data[1].matchAll(/<row\b([^>]*?)(?:\/>|>([\s\S]*?)<\/row>)/g)].map(m=>{
      const cells=[];for(const c of (m[2]||'').matchAll(/<c\b([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g)){
        const a=attrs(c[1]),ref=/^([A-Z]+)[1-9][0-9]*$/.exec(a.r||'');if(!ref)throw Error('Invalid workbook cell address.');
        const col=columnIndex(ref[1]);if(cells[col])throw Error('Duplicate workbook cell.');cells[col]={attributes:a,body:c[2]||''};
      }
      return {attributes:attrs(m[1]),cells};
    });
    if(!rows.length)throw Error('Missing workbook header.');
    const headers=rows[0].cells.map(cellValue);if(new Set(headers).size!==headers.length||headers.some(x=>!x))throw Error('Invalid workbook column headers.');
    return {sheet,headers,header:rows[0],rows:rows.slice(1)};
  }
  function tableValue(table,row,label){const index=table.headers.indexOf(label);return index<0?'':cellValue(row.cells[index]);}
  function replaceText(table,row,label,value){const index=table.headers.indexOf(label);if(index>=0)row.cells[index]=textCell(value,row.cells[index]);}
  function renderTable(table,rows,columns=table.headers.map((label,index)=>({label,index}))){
    const tagAttrs=a=>Object.entries(a).map(([k,v])=>' '+k+'="'+escapeXml(v)+'"').join('');
    const cell=(c,j,i)=>{if(!c)return '';const a={...c.attributes,r:columnName(j)+(i+1)};return '<c'+tagAttrs(a)+(c.body?'>'+c.body+'</c>':'/>');};
    const header={...table.header,cells:table.header.cells.slice()};
    const data=[header,...rows].map((r,i)=>'<row'+tagAttrs({...r.attributes,r:i+1})+'>'+columns.map((c,j)=>cell(i===0?textCell(c.label,r.cells[c.index]):r.cells[c.index],j,i)).join('')+'</row>').join('');
    let sheet=table.sheet.replace(/<sheetData\b[^>]*>[\s\S]*?<\/sheetData>/,'<sheetData>'+data+'</sheetData>');
    const ref='A1:'+columnName(Math.max(0,columns.length-1))+(rows.length+1);
    sheet=sheet.replace(/<dimension\b[^>]*\/>/,'<dimension ref="'+ref+'"/>');
    sheet=sheet.replace(/<autoFilter\b[^>]*(?:\/>|>[\s\S]*?<\/autoFilter>)/,'<autoFilter ref="'+ref+'"/>');
    const oldCols=[...table.sheet.matchAll(/<col\b([^>]*?)\/>/g)].map(m=>attrs(m[1]));
    if(oldCols.length){const cols=columns.map((c,j)=>{const source=oldCols.find(x=>Number(x.min)<=c.index+1&&Number(x.max)>=c.index+1)||{width:22,customWidth:1};return '<col'+tagAttrs({...source,min:j+1,max:j+1})+'/>';}).join('');sheet=sheet.replace(/<cols\b[^>]*>[\s\S]*?<\/cols>/,'<cols>'+cols+'</cols>');}
    return sheet;
  }
  function pairKey(table,row){return JSON.stringify([tableValue(table,row,'Comparison'),tableValue(table,row,'Role').replace(/^Neo-/,''),tableValue(table,row,'Neo form ID'),tableValue(table,row,'Spanning form ID')]);}
  function splitIDs(value){return value.split(/;\s*/).filter(Boolean);}
  function publicColumns(table,rows,{dropInternal=false,comparisonInfo=new Map()}={}){
    const excluded=new Set(dropInternal?['Control','Treatment','Control total','Treatment total','Site IDs','Export key','Passes default coverage']:[]);
    const comparisons=new Set(rows.map(r=>tableValue(table,r,'Comparison'))),name=comparisons.size===1?[...comparisons][0]:null,cp=name?comparisonInfo.get(name):null;
    return table.headers.flatMap((label,index)=>excluded.has(label)?[]:[{index,label:cp&&label==='Control measured'?cp.control+' measured (of '+cp.controlTotal+')':cp&&label==='Treatment measured'?cp.treatment+' measured (of '+cp.treatmentTotal+')':label}]);
  }
  async function filterSheets(view,load){
    if(!view||!Array.isArray(view.rows)||!view.settings||typeof view.settings!=='object')throw Error('Current view is unavailable.');
    const requested=new Map();for(const row of view.rows){if(typeof row?.key!=='string'||typeof row.pattern!=='string'||requested.has(row.key))throw Error('Invalid or duplicate current-view row.');requested.set(row.key,row);}
    const source=dataTable(await load('__Results')),byKey=new Map();
    for(const row of source.rows){const key=tableValue(source,row,'Export key');if(!key||byKey.has(key))throw Error('Invalid canonical result key.');byKey.set(key,row);}
    const picked=[...requested].map(([key,v])=>{const original=byKey.get(key);if(!original)throw Error('Current view does not match the saved workbook.');const row={...original,cells:original.cells.slice()};replaceText(source,row,'Pattern',v.pattern);return row;});
    const mappings=new Set(),forms=new Set(),sites=new Set(),comparisonInfo=new Map();
    for(const row of picked){const comp=tableValue(source,row,'Comparison'),role=tableValue(source,row,'Role');
      for(const site of splitIDs(tableValue(source,row,'Site IDs'))){sites.add(site);mappings.add(JSON.stringify([comp,role,site]));}
      for(const column of ['Neo form ID','Spanning form ID'])forms.add(tableValue(source,row,column));
      const info={control:tableValue(source,row,'Control'),treatment:tableValue(source,row,'Treatment'),controlTotal:tableValue(source,row,'Control total'),treatmentTotal:tableValue(source,row,'Treatment total')};
      if(comparisonInfo.has(comp)&&JSON.stringify(comparisonInfo.get(comp))!==JSON.stringify(info))throw Error('Inconsistent comparison metadata.');comparisonInfo.set(comp,info);
    }
    const patterns=new Map((view.pairPatterns||[]).map(x=>[x.key,x.pattern]));
    const all=dataTable(await load('All pairs')),allRows=all.rows.filter(row=>splitIDs(tableValue(all,row,'Site IDs')).some(site=>mappings.has(JSON.stringify([tableValue(all,row,'Comparison'),tableValue(all,row,'Role'),site]))));
    for(const row of allRows){for(const column of ['Neo form ID','Spanning form ID'])forms.add(tableValue(all,row,column));const pattern=patterns.get(pairKey(all,row));if(pattern!==undefined)replaceText(all,row,'Pattern',pattern);}
    const transformed=new Map([['Results',renderTable(source,picked,publicColumns(source,picked,{dropInternal:true,comparisonInfo}))],['All pairs',renderTable(all,allRows,publicColumns(all,allRows,{comparisonInfo}))]]);
    return async function transform(name,sheet){
      if(transformed.has(name))return transformed.get(name);
      const table=dataTable(sheet);
      if(name==='Peptides'||name==='Other candidates')return renderTable(table,table.rows.filter(row=>forms.has(tableValue(table,row,'Form ID'))));
      if(name==='Annotations')return renderTable(table,table.rows.filter(row=>sites.has(tableValue(table,row,'site_id'))));
      if(name==='Settings'){
        const rows=table.rows.filter(row=>tableValue(table,row,'Section')!=='Export');
        const labels={comparison:'Comparison',control:'Control',treatment:'Treatment',controlTotal:'Control total',treatmentTotal:'Treatment total',threshold:'Ratio threshold',coverage:'Display coverage k',analysisSelectionK:'Frozen selection k',search:'Search',role:'Role filter',support:'Terminal evidence filter',spanning:'Spanning filter',referenceOnly:'Reference filter',protease:'Protease filter',proteinRecord:'Protein record filter',pattern:'Pattern filter',sortField:'Sort field',sortDirection:'Sort direction',rowCount:'Result rows',scope:'Scope'};
        const snapshot={...view.settings,scope:'Current view',rowCount:picked.length};
        for(const [key,value] of Object.entries(snapshot)){
          const row={attributes:{},cells:[]};for(const [label,text] of [['Section','Export'],['Setting',labels[key]||key],['Value',value==null?'':typeof value==='object'?JSON.stringify(value):String(value)]]){const i=table.headers.indexOf(label);if(i<0)throw Error('Invalid Settings table.');row.cells[i]=textCell(text,table.rows[0]?.cells[i]);}rows.push(row);
        }
        const row={attributes:{},cells:[]};for(const [label,text] of [['Section','Export'],['Setting','Related tables scope'],['Value','Peptides includes forms from all alternative pairs at the retained comparison/role/sites. Annotations includes retained sites. Sample quantities retain their original columns. Tests and BH values are unchanged.']]){const i=table.headers.indexOf(label);row.cells[i]=textCell(text,table.rows[0]?.cells[i]);}rows.push(row);
        return renderTable(table,rows);
      }
      return sheet;
    };
  }
  async function select(workbook,selected,options={}){
    if(!workbook||!Array.isArray(workbook.sheets)||!workbook.sheets.length||new Set(workbook.sheets).size!==workbook.sheets.length||
      workbook.sheets.some(name=>typeof name!=='string'||!name))throw Error('Invalid workbook table list.');
    if(!Array.isArray(selected)||!selected.length||new Set(selected).size!==selected.length||selected.some(name=>!workbook.sheets.includes(name)))throw Error('Select at least one available table, without duplicates.');
    if(options.view&&!selected.includes('Settings'))selected=[...selected,'Settings'];
    const zip=archive(base64(workbook.base64)),book=xmlText(await zip.read('xl/workbook.xml')),
      rels=relationships(xmlText(await zip.read('xl/_rels/workbook.xml.rels')),'xl/workbook.xml'),sheets=new Map();
    for(const m of book.matchAll(/<sheet\b([^>]*?)\/?\s*>/g)){const a=attrs(m[1]),rel=rels.get(a['r:id']);
      if(!a.name||sheets.has(a.name)||!rel||rel.type!==REL+'/worksheet')throw Error('Invalid workbook worksheet relationship.');sheets.set(a.name,rel.path);}
    const internal=workbook.internal_sheets||[],declared=[...workbook.sheets,...internal];
    if(!Array.isArray(internal)||internal.some(n=>n!=='__Results')||new Set(declared).size!==declared.length||declared.length!==sheets.size||declared.some(name=>!sheets.has(name)))throw Error('Workbook table list does not match its contents.');
    let strings=[];const shared=[...rels.values()].find(r=>r.type===REL+'/sharedStrings');
    if(shared){const xml=xmlText(await zip.read(shared.path));strings=[...xml.matchAll(/<si\b[^>]*?(?:\/>|>([\s\S]*?)<\/si>)/g)].map(m=>m[1]??'');}
    const loaded=new Map();
    async function load(name){
      if(loaded.has(name))return loaded.get(name);
      if(!sheets.has(name))throw Error('Filtered export requires a refreshed report with '+name+'.');
      let sheet=xmlText(await zip.read(sheets.get(name)));
      if(/<(?:f|drawing|legacyDrawing|hyperlink|tablePart|oleObject|control)\b/.test(sheet))throw Error('Choose tables supports ordinary Atelente data worksheets only.');
      sheet=sheet.replace(/<c\b([^>]*?)(?:\/>|>([\s\S]*?)<\/c>)/g,(whole,attributes,contents)=>{
        if(contents===undefined||attrs(attributes).t!=='s')return whole;
        const m=contents.match(/^\s*<v>\s*(\d+)\s*<\/v>\s*$/);if(!m||Number(m[1])>=strings.length)throw Error('Invalid workbook shared-string index.');
        return '<c'+attributes.replace(/\bt\s*=\s*(["'])s\1/,'t="inlineStr"')+'><is>'+strings[Number(m[1])]+'</is></c>';
      });loaded.set(name,sheet);return sheet;
    }
    const transform=options.view?await filterSheets(options.view,load):null;
    const declaration='<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',parts=[],sheetTags=[],sheetRels=[],overrides=[],defined=[];
    for(const [index,name] of selected.entries()){
      let sheet=await load(name);if(transform)sheet=await transform(name,sheet);
      // Printer-specific binary relationships are not part of data export.
      sheet=sheet.replace(/<pageSetup\b[^>]*>/g,tag=>tag.replace(/\s+r:id\s*=\s*(["']).*?\1/g,''));
      if(/\br:id\s*=/.test(sheet))throw Error('Unsupported worksheet relationship in selected table.');
      sheet=sheet.replace(/\btabSelected\s*=\s*(["']).*?\1/g,'tabSelected="'+(index===0?'true':'false')+'"');
      const number=index+1,target='worksheets/sheet'+number+'.xml';parts.push(['xl/'+target,sheet]);
      sheetTags.push('<sheet name="'+escapeXml(name)+'" sheetId="'+number+'" r:id="rId'+number+'"/>');
      sheetRels.push('<Relationship Id="rId'+number+'" Type="'+REL+'/worksheet" Target="'+target+'"/>');
      overrides.push('<Override PartName="/xl/'+target+'" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>');
      const filter=sheet.match(/<autoFilter\b([^>]*?)\/?\s*>/);if(filter){const ref=attrs(filter[1]).ref;
        if(!/^[A-Z]+[0-9]+(?::[A-Z]+[0-9]+)?$/.test(ref))throw Error('Unsupported worksheet filter reference.');
        const absolute=ref.replace(/([A-Z]+)([0-9]+)/g,'$$$1$$$2');
        defined.push('<definedName name="_xlnm._FilterDatabase" localSheetId="'+index+'" hidden="1">'+escapeXml("'"+name.replace(/'/g,"''")+"'!"+absolute)+'</definedName>');}
    }
    for(const type of ['styles','theme']){const source=[...rels.values()].find(r=>r.type===REL+'/'+type);if(!source)continue;
      const target=type==='styles'?'styles.xml':'theme/theme1.xml';parts.push(['xl/'+target,xmlText(await zip.read(source.path))]);
      sheetRels.push('<Relationship Id="rId'+type+'" Type="'+REL+'/'+type+'" Target="'+target+'"/>');
      overrides.push('<Override PartName="/xl/'+target+'" ContentType="'+(type==='styles'?'application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml':'application/vnd.openxmlformats-officedocument.theme+xml')+'"/>');}
    const property=book.match(/<workbookPr\b([^>]*?)\/?\s*>/),date1904=property&&['1','true'].includes(attrs(property[1]).date1904);
    parts.push(['xl/workbook.xml',declaration+'<workbook xmlns="'+MAIN+'" xmlns:r="'+REL+'"><workbookPr date1904="'+(date1904?'true':'false')+'"/><bookViews><workbookView activeTab="0"/></bookViews><sheets>'+sheetTags.join('')+'</sheets>'+(defined.length?'<definedNames>'+defined.join('')+'</definedNames>':'')+'</workbook>']);
    parts.push(['xl/_rels/workbook.xml.rels',declaration+'<Relationships xmlns="'+PKG+'">'+sheetRels.join('')+'</Relationships>']);
    parts.push(['_rels/.rels',declaration+'<Relationships xmlns="'+PKG+'"><Relationship Id="rId1" Type="'+REL+'/officeDocument" Target="xl/workbook.xml"/></Relationships>']);
    parts.push(['[Content_Types].xml',declaration+'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'+overrides.join('')+'</Types>']);
    return storeZip(parts);
  }
  function mount(D,document,options={}){
    const ids=Object.fromEntries(['downloads','tableChooser','tableChoices','downloadChosen','downloadStatus'].map(id=>[id,document.getElementById(id)]));
    if(Object.values(ids).some(x=>!x))throw Error('Missing workbook download controls.');
    const create=(tag,text,parent)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;parent.append(e);return e;};
    ids.downloads.replaceChildren();ids.tableChoices.replaceChildren();ids.tableChooser.hidden=true;ids.downloadStatus.textContent='';ids.downloadChosen.onclick=null;
    const link=(label,href)=>{if(typeof href!=='string'||!href)return;safePath(href);const a=create('a',label,ids.downloads);a.href=href;a.download=href.split('/').at(-1);return a;};
    link('Download results — Excel',D.excel_file);
    const choose=create('button','Choose tables…',ids.downloads);choose.type='button';choose.setAttribute('aria-expanded','false');choose.setAttribute('aria-controls','tableChooser');
    link('Download full audit — ZIP',D.audit_file);
    const workbook=D.export_workbook,available=Array.isArray(workbook?.sheets)&&typeof workbook?.base64==='string',checks=[];
    choose.disabled=!available;ids.downloadChosen.type='button';ids.downloadChosen.textContent='Download selected tables';ids.downloadChosen.disabled=!available;
    if(!available){choose.title='Custom tables are unavailable in this report.';return {downloadCurrent:async()=>{throw Error('Excel export is unavailable in this report.');}};}
    const canFilter=typeof options.getView==='function'&&workbook.internal_sheets?.includes('__Results');
    let scope=null;if(canFilter){const label=create('label','Export scope ',ids.tableChoices);scope=create('select',undefined,label);scope.id='exportScope';for(const [value,text] of [['all','All results'],['current','Current view']]){const option=create('option',text,scope);option.value=value;}scope.value='all';}
    for(const name of workbook.sheets){const label=create('label',undefined,ids.tableChoices),check=create('input',undefined,label);check.type='checkbox';check.value=name;
      check.checked=['Results','Peptides','Settings'].includes(name);create('span',name,label);checks.push(check);
      check.onchange=()=>{ids.downloadChosen.disabled=!checks.some(c=>c.checked);ids.downloadStatus.textContent='';};}
    ids.downloadChosen.disabled=!checks.some(c=>c.checked);
    if(scope)scope.onchange=()=>{
      const settings=checks.find(c=>c.value==='Settings');
      if(settings){settings.disabled=scope.value==='current';if(settings.disabled)settings.checked=true;}
      ids.downloadChosen.disabled=!checks.some(c=>c.checked);ids.downloadStatus.textContent='';
    };
    choose.onclick=()=>{ids.tableChooser.hidden=!ids.tableChooser.hidden;choose.setAttribute('aria-expanded',String(!ids.tableChooser.hidden));};
    async function download(selected,current){
      ids.downloadChosen.disabled=true;ids.downloadStatus.textContent='Preparing Excel…';
      try{
        // Capture UI state once at click time; later UI changes cannot alter it.
        const view=current?options.getView():null;
        if(current&&!canFilter)throw Error('Current-view export requires a refreshed report.');
        const bytes=await select(workbook,selected,view?{view}:{}),url=URL.createObjectURL(new Blob([bytes],{type:'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'}));
        const anchor=document.createElement('a');anchor.href=url;anchor.download=current?'atelente_current_view.xlsx':'atelente_selected_tables.xlsx';anchor.click();setTimeout(()=>URL.revokeObjectURL(url),30000);
        ids.downloadStatus.textContent='Excel prepared · '+(view?view.rows.length+' results · ':'')+(new Set([...selected,...(view?['Settings']:[])])).size+' tables';
      }
      catch(error){ids.downloadStatus.textContent='Could not create Excel: '+error.message;}
      finally{ids.downloadChosen.disabled=!checks.some(c=>c.checked);}
    }
    ids.downloadChosen.onclick=()=>download(checks.filter(c=>c.checked).map(c=>c.value),scope?.value==='current');
    return Object.freeze({downloadCurrent:()=>download(['Results','Peptides','Settings'].filter(s=>workbook.sheets.includes(s)),true)});
  }
  return Object.freeze({select,mount});
});
