/* Versioned offline report transport, including lossless column dictionaries. */
(function(root,factory){
  if(typeof module==='object'&&module.exports)module.exports=factory(require('node:zlib'));
  else root.ReportPayload=factory(null);
})(typeof globalThis!=='undefined'?globalThis:this,function(zlib){
  'use strict';
  const SCHEMA='atelente-report-envelope-1',THRESHOLD=1024*1024;
  const DATA_SCHEMA='atelente-report-data-2';
  function unpack(value){
    if(!value||typeof value!=='object')return value;
    if(value.__atelente_table__===1){
      const n=value.nrow,columns=value.columns;
      if(!Number.isSafeInteger(n)||n<0||!columns||typeof columns!=='object')throw Error('Invalid report table.');
      const rows=Array.from({length:n},()=>({}));
      for(const [name,column] of Object.entries(columns)){
        if(!['integer','numeric','logical','character'].includes(column.type)||!Array.isArray(column.values)||column.values.length!==n)throw Error('Invalid report column.');
        const dict=column.dictionary;
        if(dict!==undefined&&!Array.isArray(dict))throw Error('Invalid report dictionary.');
        for(let i=0;i<n;i++){
          let v=column.values[i];
          if(dict&&v!==null){if(!Number.isInteger(v)||v<0||v>=dict.length)throw Error('Invalid report dictionary index.');v=dict[v];}
          if(v!==null&&!(column.type==='character'?typeof v==='string':column.type==='logical'?typeof v==='boolean':column.type==='integer'?Number.isInteger(v):typeof v==='number'&&Number.isFinite(v)))throw Error('Invalid report column value.');
          Object.defineProperty(rows[i],name,{value:v,writable:true,enumerable:true,configurable:true});
        }
      }
      return rows;
    }
    if(Array.isArray(value))return value.map(unpack);
    for(const name of Object.keys(value))value[name]=unpack(value[name]);
    return value;
  }
  function validate(value){
    if(value?.payload_schema!==undefined&&value.payload_schema!==DATA_SCHEMA)throw Error('Unsupported report data schema. Use a compatible Atelente version.');
    return value?.payload_schema===DATA_SCHEMA?unpack(value):value;
  }
  function envelope(value){
    if(!value||value.schema!==SCHEMA)return false;
    if(value.encoding!=='deflate-base64'||!Number.isSafeInteger(value.uncompressed_bytes)||value.uncompressed_bytes<0||
      typeof value.data!=='string'||value.data.length%4!==0||/[^A-Za-z0-9+/=]/.test(value.data)||/=/.test(value.data.slice(0,-2)))
      throw Error('Invalid compressed report payload.');
    return true;
  }
  function decoded(bytes,value){
    if(bytes.byteLength!==value.uncompressed_bytes)throw Error('Report payload length check failed.');
    return validate(JSON.parse(new TextDecoder('utf-8',{fatal:true}).decode(bytes)));
  }
  function parseSync(text){
    const value=JSON.parse(text);
    if(!envelope(value))return validate(value);
    if(!zlib)throw Error('Synchronous report decoding is available only in Node.');
    return decoded(zlib.inflateSync(Buffer.from(value.data,'base64')),value);
  }
  async function parse(text){
    const value=JSON.parse(text);
    if(!envelope(value))return validate(value);
    if(typeof DecompressionStream!=='function')throw Error('This compressed report requires a browser with DecompressionStream support. Open it in a current Chrome, Edge, Firefox or Safari version.');
    const bytes=Uint8Array.from(atob(value.data),c=>c.charCodeAt(0));
    const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));
    return decoded(new Uint8Array(await new Response(stream).arrayBuffer()),value);
  }
  function encodeSync(json,threshold=THRESHOLD){
    if(!zlib)throw Error('Synchronous report encoding is available only in Node.');
    const bytes=Buffer.from(json,'utf8');
    if(bytes.length<threshold)return json;
    const text=JSON.stringify({schema:SCHEMA,encoding:'deflate-base64',uncompressed_bytes:bytes.length,data:zlib.deflateSync(bytes).toString('base64')});
    return Buffer.byteLength(text)<bytes.length?text:json;
  }
  return {SCHEMA,DATA_SCHEMA,THRESHOLD,parse,parseSync,encodeSync};
});
