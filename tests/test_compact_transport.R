script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/')
source(file.path(root,'R/report.R'),encoding='UTF-8')
x<-list(payload_schema='atelente-report-data-2', rows=data.frame(id=c('001','<script>','T\u00fcrk\u00e7e',NA_character_),value=c(pi,0,NA_real_,1e-100),flag=c(TRUE,FALSE,NA,TRUE),count=c(1L,2L,NA_integer_,4L)), empty=data.frame(id=character(),value=numeric()),single=data.frame(value=pi), repeated=data.frame(id=rep(c('\u03b1','001',NA_character_),30),v=seq_len(90L)))
packed<-as.character(jsonlite::toJSON(pm_report_pack_tables(x),auto_unbox=TRUE,digits=I(17),na='null',null='null'))
y<-pm_report_decode_json(packed)
stopifnot(identical(x,y))
stopifnot(inherits(try(pm_report_decode_json('{"payload_schema":"future"}'),silent=TRUE),'try-error'))
# Lists of tables must not be mistaken for one ordinary JSON row table. Likewise,
# a character dictionary contains numeric codes even with simplifyVector=FALSE.
x$nested <- list(tables=list(data.frame(a=1L),data.frame(a=2L)),metadata=list(list(name='first'),list(name='second')))
x$zero_columns <- structure(list(),class='data.frame',row.names=1:2)
x$all_missing <- data.frame(text=rep(NA_character_,40L),value=rep(NA_real_,40L))
x$precision <- data.frame(value=c(1+.Machine$double.eps,1-.Machine$double.eps/2,
  exp(1),.Machine$double.xmin,.Machine$double.xmin/2,.Machine$double.xmax))
x$names <- data.frame(matrix(c('literal','001'),nrow=1L),check.names=FALSE)
names(x$names) <- c('__proto__','constructor')
packed <- as.character(jsonlite::toJSON(pm_report_pack_tables(x),auto_unbox=TRUE,digits=I(17),na='null',null='null'))
stopifnot(identical(x,pm_report_decode_json(packed)),identical(x,pm_report_decode_json(packed,simplifyVector=FALSE)))
envelope <- pm_report_encode_json(packed,threshold=0)
stopifnot(identical(x,pm_report_decode_json(envelope)))
malformed <- pm_report_pack_tables(x)
malformed$rows$columns$value$values <- I(1)
stopifnot(inherits(try(pm_report_decode_json(jsonlite::toJSON(malformed,auto_unbox=TRUE)),silent=TRUE),'try-error'))
malformed <- pm_report_pack_tables(x)
malformed$repeated$columns$id$values[1] <- 999L
stopifnot(inherits(try(pm_report_decode_json(jsonlite::toJSON(malformed,auto_unbox=TRUE)),silent=TRUE),'try-error'))
malformed <- pm_report_pack_tables(x); malformed$rows$columns$value$values <- I(rep('1',4L))
stopifnot(inherits(try(pm_report_decode_json(jsonlite::toJSON(malformed,auto_unbox=TRUE)),silent=TRUE),'try-error'))
legacy <- '{"rows":[{"id":"001","value":0},{"id":"002","value":null}],"labels":["A","B"]}'
stopifnot(identical(pm_report_decode_json(legacy),jsonlite::fromJSON(legacy)),
  identical(pm_report_decode_json(legacy,simplifyVector=FALSE),jsonlite::fromJSON(legacy,simplifyVector=FALSE)))
# Browser decoding must preserve the exact same numeric values and literal IDs.
# This is a codec test, not a claim of real-browser rendering/memory validation.
node <- Sys.which('node')
if(nzchar(node)) {
  folder <- tempfile('atelente-codec-');dir.create(folder)
  writeBin(charToRaw(enc2utf8(packed)),file.path(folder,'packed.json'))
  writeBin(charToRaw(enc2utf8(envelope)),file.path(folder,'envelope.json'))
  writeBin(charToRaw(enc2utf8(as.character(jsonlite::toJSON(x,auto_unbox=TRUE,dataframe='rows',digits=I(17),na='null',null='null')))),file.path(folder,'expected.json'))
  script_path <- file.path(folder,'check.cjs')
  writeLines(c("const fs=require('node:fs'),assert=require('node:assert/strict'),path=require('node:path');",
    "const codec=require(process.argv[2]),folder=process.argv[3],read=n=>fs.readFileSync(path.join(folder,n),'utf8');",
    "const expected=JSON.parse(read('expected.json'));",
    "assert.deepStrictEqual(codec.parseSync(read('packed.json')),expected);assert.deepStrictEqual(codec.parseSync(read('envelope.json')),expected);",
    "assert.equal(codec.parseSync(read('packed.json')).rows[0].value,Math.PI);",
    "assert.equal(codec.parseSync(read('packed.json')).rows[1].value,0);assert.equal(codec.parseSync(read('packed.json')).rows[2].value,null);",
    "assert.equal(codec.parseSync(read('packed.json')).precision[0].value,1+Number.EPSILON);",
    "assert.throws(()=>codec.parseSync('{\"payload_schema\":\"future\"}'),/schema/);",
    "const malformed=JSON.parse(read('packed.json'));malformed.repeated.columns.id.values[0]=999;assert.throws(()=>codec.parseSync(JSON.stringify(malformed)),/dictionary/);",
    "const wrongType=JSON.parse(read('packed.json'));wrongType.rows.columns.value.values[0]='1';assert.throws(()=>codec.parseSync(JSON.stringify(wrongType)),/column value/);",
    "console.log('R-to-JS exact values, typed nested tables, legacy envelope and malformed dictionaries passed.');"),script_path,useBytes=TRUE)
  result <- system2(node,c(shQuote(script_path),shQuote(file.path(root,'inst','templates','report_payload.js')),shQuote(folder)))
  stopifnot(result==0L)
}
cat('Column transport R roundtrip, null/zero, literal IDs, nested lists and schema guards passed.\n')
