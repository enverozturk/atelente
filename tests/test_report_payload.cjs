// Lossless transport, browser decoder, backward compatibility and visible failure contract.
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm');
const C=require('../inst/templates/report_payload.js');
const input={unicode:'Türkçe α β',literal:'</script>&',zero:0,missing:null,numeric:Math.PI,
  rows:Array.from({length:30000},(_,i)=>({id:i,form:'PEPTIDE',value:i/7,missing:null}))};
const json=JSON.stringify(input),packed=C.encodeSync(json);
assert.equal(JSON.parse(packed).schema,C.SCHEMA);
assert.ok(Buffer.byteLength(packed)<Buffer.byteLength(json)/3);
assert.deepEqual(C.parseSync(packed),input);
assert.deepEqual(C.parseSync(json),input);
assert.equal(C.encodeSync('{"small":true}'),'{"small":true}');
const bad=JSON.parse(packed);bad.uncompressed_bytes++;
assert.throws(()=>C.parseSync(JSON.stringify(bad)),/length check/);
bad.encoding='unknown';assert.throws(()=>C.parseSync(JSON.stringify(bad)),/Invalid/);
const code=fs.readFileSync(require.resolve('../inst/templates/report_payload.js'),'utf8');
assert.ok(!/<\/script\s*>/i.test(code));
(async()=>{
  // Actual browser branch: no Node require/module provided.
  const sandbox={DecompressionStream,TextDecoder,Blob,Response,Uint8Array,atob};
  vm.runInNewContext(code,sandbox);
  assert.equal(JSON.stringify(await sandbox.ReportPayload.parse(packed)),json);
  assert.equal(JSON.stringify(await sandbox.ReportPayload.parse(json)),json);
  const unsupported={};vm.runInNewContext(code,unsupported);
  await assert.rejects(()=>unsupported.ReportPayload.parse(packed),/DecompressionStream support/);
  const template=fs.readFileSync(require.resolve('../inst/templates/cleavage_report_template.html'),'utf8');
  assert.ok(template.includes("notice.setAttribute('role','alert')"));
  for(const file of process.argv.slice(2)){
    const html=fs.readFileSync(file,'utf8'),text=html.match(/<script type="application\/json" id="data">([\s\S]*?)<\/script>/)[1];
    const canonical=JSON.stringify(C.parseSync(text)),repacked=C.encodeSync(canonical);
    assert.equal(JSON.stringify(C.parseSync(repacked)),canonical);
    console.log(JSON.stringify({file,input_bytes:Buffer.byteLength(text),packed_bytes:Buffer.byteLength(repacked)}));
  }
  console.log('Lossless report payload, native browser decode and legacy compatibility passed.');
})().catch(error=>{console.error(error);process.exitCode=1});
