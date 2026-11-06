/* Presentation summaries only. No model fitting, p-value selection or imputation. */
(function(root,factory){
  const api=factory();
  if(typeof module==='object'&&module.exports)module.exports=api;
  else root.PairSummary=api;
})(typeof globalThis!=='undefined'?globalThis:this,function(){
  'use strict';
  const SELECTION_POLICY='ends_coverage_overlap_intensity_v3';
  const finite=Number.isFinite;
  const mean=values=>values.length?values.reduce((sum,x)=>sum+x,0)/values.length:null;
  const difference=(a,b)=>finite(a)&&finite(b)?a-b:null;
  const positive=(r,role)=>r?.[role+'_status']==='observed'&&finite(r[role+'_value'])&&r[role+'_value']>0;
  const close=(a,b)=>finite(a)&&finite(b)?Math.abs(a-b)<=1e-10*Math.max(1,Math.abs(a),Math.abs(b)):a===b;
  const literal=(a,b)=>String(a)<String(b)?-1:String(a)>String(b)?1:0;
  const pairKey=r=>JSON.stringify([r.site_id,r.role,r.neo_peptide_id,r.spanning_peptide_id]);
  const measurementKey=r=>JSON.stringify([r.comparison,r.role,r.neo_peptide_id,r.spanning_peptide_id]);
  function selectedSamples(samples,comparison){
    const chosen=(samples||[]).filter(s=>s.Use==='Yes'&&(!comparison||[comparison.Control,comparison.Treatment].includes(s.Condition)));
    if(new Set(chosen.map(s=>s['Peptide column'])).size!==chosen.length)throw new Error('Duplicate included sample columns.');
    return chosen;
  }
  function defaultK(samples,comparison){
    if(!comparison)return 0;
    const chosen=selectedSamples(samples,comparison);
    const n=[comparison.Control,comparison.Treatment].map(c=>chosen.filter(s=>s.Condition===c).length);
    return Math.min(3,...n);
  }
  function conditionSummary(rows,samples,diagnostics){
    const lookup=new Map(rows.map(r=>[r.source_column,r]));
    const members=samples.map(s=>lookup.get(s['Peptide column']));
    const neo=members.filter(r=>positive(r,'neo')),span=members.filter(r=>positive(r,'spanning'));
    const co=members.filter(r=>positive(r,'neo')&&positive(r,'spanning')&&finite(r.log2_ratio));
    for(const row of members){
      if(!row){diagnostics.push('missing_sample_row');continue;}
      const both=positive(row,'neo')&&positive(row,'spanning');
      if(both&&!close(row.log2_ratio,Math.log2(row.neo_value)-Math.log2(row.spanning_value)))diagnostics.push('sample_ratio_mismatch');
      if(!both&&finite(row.log2_ratio))diagnostics.push('ratio_without_positive_components');
    }
    return {total:members.length,co:co.length,neo:neo.length,spanning:span.length,
      neoMissing:members.filter(r=>r?.neo_status==='missing').length,
      spanningMissing:members.filter(r=>r?.spanning_status==='missing').length,
      neoZero:members.filter(r=>r?.neo_status==='zero').length,
      spanningZero:members.filter(r=>r?.spanning_status==='zero').length,
      meanNeoLog2:mean(co.map(r=>Math.log2(r.neo_value))),
      meanSpanLog2:mean(co.map(r=>Math.log2(r.spanning_value))),
      meanLog2Ratio:mean(co.map(r=>r.log2_ratio)),
      meanAllNeoLog2:mean(neo.map(r=>Math.log2(r.neo_value))),
      meanAllSpanLog2:mean(span.map(r=>Math.log2(r.spanning_value)))};
  }
  function summarizePairs(rows,samples,comparison,contrasts=[]){
    const grouped=new Map(),exports=new Map();
    for(const row of rows||[]){
      const key=pairKey(row);
      if(!grouped.has(key))grouped.set(key,[]);
      grouped.get(key).push(row);
    }
    for(const row of contrasts||[])if(comparison&&row.comparison===comparison.Name){
      const key=pairKey(row);
      if(exports.has(key))throw new Error('Duplicate exported pair contrast.');
      exports.set(key,row);
    }
    const selected=selectedSamples(samples,comparison);
    return [...grouped].map(([key,observations])=>{
      if(new Set(observations.map(r=>r.source_column)).size!==observations.length)throw new Error('Duplicate pair/sample ratio rows.');
      const r=observations[0],diagnostics=[];
      const control=conditionSummary(observations,comparison?selected.filter(s=>s.Condition===comparison.Control):[],diagnostics);
      const treatment=conditionSummary(observations,comparison?selected.filter(s=>s.Condition===comparison.Treatment):[],diagnostics);
      let effect=difference(treatment.meanLog2Ratio,control.meanLog2Ratio);
      const exported=exports.get(key);
      if(exported){
        if(!close(effect,exported.log2_fc??null))diagnostics.push('exported_effect_mismatch');
        for(const [side,counts] of [['control',control],['treatment',treatment]]){
          for(const [suffix,field] of [['samples','total'],['neo','neo'],['spanning','spanning'],['usable_pairs','co']]){
            if(exported[side+'_'+suffix]!==undefined&&exported[side+'_'+suffix]!==counts[field])diagnostics.push('exported_coverage_mismatch');
          }
        }
        effect=exported.log2_fc??null;
      }
      const neoFC=difference(treatment.meanNeoLog2,control.meanNeoLog2);
      const spanFC=difference(treatment.meanSpanLog2,control.meanSpanLog2);
      if(finite(effect)&&!close(effect,neoFC-spanFC))diagnostics.push('component_effect_mismatch');
      const valid=!diagnostics.length;
      const computedFold=finite(effect)?Math.pow(2,effect):null;
      const fold=exported?(exported.fold_change??null):finite(computedFold)&&computedFold>0?computedFold:null;
      if(exported&&fold!==null&&finite(effect)&&!close(Math.log2(fold),effect))diagnostics.push('exported_fold_mismatch');
      const arithmeticValid=!diagnostics.length;
      // Pool included common samples without condition weighting. Stable column
      // order makes floating-point accumulation independent of input row order.
      const selectedColumns=new Set(selected.map(s=>s['Peptide column']));
      const common=observations.filter(row=>selectedColumns.has(row.source_column)&&positive(row,'neo')&&positive(row,'spanning')&&finite(row.log2_ratio))
        .sort((a,b)=>literal(a.source_column,b.source_column));
      const commonLogIntensity=arithmeticValid?mean(common.map(row=>(Math.log2(row.neo_value)+Math.log2(row.spanning_value))/2)):null;
      return {key,site_id:r.site_id,record_id:r.record_id,accession:r.accession,bond_left:r.bond_left,bond_right:r.bond_right,
        role:r.role,neo_peptide_id:r.neo_peptide_id,spanning_peptide_id:r.spanning_peptide_id,
        comparison:comparison?.Name??null,control_name:comparison?.Control??null,treatment_name:comparison?.Treatment??null,
        control,treatment,neo_log2_fc:arithmeticValid?neoFC:null,spanning_log2_fc:arithmeticValid?spanFC:null,
        log2_fc:arithmeticValid?effect:null,fold_change:arithmeticValid?fold:null,
        neo_all_log2_fc:valid?difference(treatment.meanAllNeoLog2,control.meanAllNeoLog2):null,
        spanning_all_log2_fc:valid?difference(treatment.meanAllSpanLog2,control.meanAllSpanLog2):null,
        selection_k:comparison?Math.max(1,defaultK(samples,comparison)):0,common_log2_intensity:commonLogIntensity,
        exported_status:exported?.status??null,arithmetic_valid:arithmeticValid,diagnostics:[...new Set(diagnostics)],rows:observations};
    });
  }
  // R supplies arithmetic and ordering; the browser changes display filters only.
  function preparedPairs(rows,comparison,digestionLookup){
    return (rows||[]).filter(r=>(r.comparison??null)===(comparison?.Name??null)).map(r=>{
      const condition=side=>{
        const result={};
        for(const name of ['total','neo','spanning','neoMissing','spanningMissing','neoZero','spanningZero','meanNeoLog2','meanSpanLog2','meanLog2Ratio','meanAllNeoLog2','meanAllSpanLog2'])result[name]=r[side+'_'+name];
        result.co=r[side+'_common'];return result;
      };
      const p={key:pairKey(r),site_id:r.site_id,record_id:r.record_id,accession:r.accession,bond_left:r.bond_left,bond_right:r.bond_right,
        role:r.role,neo_peptide_id:r.neo_peptide_id,spanning_peptide_id:r.spanning_peptide_id,comparison:r.comparison,
        neo_log2_fc:r.neo_log2_fc,spanning_log2_fc:r.spanning_log2_fc,log2_fc:r.log2_fc,fold_change:r.fold_change,
        neo_all_log2_fc:r.neo_all_log2_fc,spanning_all_log2_fc:r.spanning_all_log2_fc,arithmetic_valid:r.arithmetic_valid,
        selection_k:r.selection_k,control_name:r.control,treatment_name:r.treatment,control:condition('control'),treatment:condition('treatment'),
        common_log2_intensity:r.intensity,exported_status:r.status,diagnostics:[],rows:[],
        selection:{rule_version:SELECTION_POLICY,selection_k:r.selection_k,eligible:r.eligible,
          preferred_digestion:r.preferred,other_candidate_bonds:r.overlap,common_log2_intensity:r.intensity,
          fallback:r.selection_fallback,rank:r.selection_rank}};
      Object.defineProperty(p.selection,'digestion_reasons',{enumerable:true,get(){return digestionLookup?digestionPreference(digestionLookup(p)).reasons:[];}});
      return p;
    });
  }
  function preparedOrder(pairs){
    const ranks=new Set();
    for(const p of pairs){const rank=p.selection?.rank;if(!Number.isInteger(rank)||rank<1||ranks.has(rank))throw Error('Invalid prepared representative order.');ranks.add(rank);}
    return [...pairs].sort((a,b)=>a.selection.rank-b.selection.rank);
  }
  const coverageOrder=(a,b)=>Math.min(b.control.co,b.treatment.co)-Math.min(a.control.co,a.treatment.co)||
    (b.control.co+b.treatment.co)-(a.control.co+a.treatment.co);
  const identityOrder=(a,b)=>literal(a.neo_peptide_id,b.neo_peptide_id)||
    literal(a.spanning_peptide_id,b.spanning_peptide_id)||literal(a.key,b.key);
  function legacyRankPairs(pairs){
    return [...pairs].sort((a,b)=>coverageOrder(a,b)||identityOrder(a,b));
  }
  function digestionPreference(metadata){
    const d=metadata||{},reasons=[];
    for(const name of ['neo_other_end_compatible','spanning_n_compatible','spanning_c_compatible']){
      if(d[name]!==true)reasons.push(name+(d[name]===false?'_false':'_unknown'));
    }
    // Internal missed-cleavage counts are descriptive metadata, not a ranking
    // criterion. Unknown MC alone cannot override known-compatible ends.
    return {preferred:!reasons.length,reasons};
  }
  function buildDigestionLookup(payload){
    const cleavage=payload?.cleavage||{},profile=cleavage.profile||{};
    const knownProfile=profile.enzyme==='trypsin'&&typeof profile.lysine_blocked==='boolean';
    const peptides=new Map((payload?.peptides||[]).map(p=>[p.peptide_id,p]));
    const evidence=new Map(),termini=new Map(),sitesByRecord=new Map(),sitesById=new Map(),overlapCache=new Map();
    for(const site of cleavage.sites||[]){
      sitesById.set(site.site_id,site);
      if(!sitesByRecord.has(site.record_id))sitesByRecord.set(site.record_id,new Set());
      sitesByRecord.get(site.record_id).add(site.bond_left);
    }
    for(const [record,positions] of sitesByRecord)sitesByRecord.set(record,[...positions].sort((a,b)=>a-b));
    const evidenceKey=(site,record,role,peptide)=>JSON.stringify([site,record,role,peptide]);
    const terminalKey=(placement,role)=>JSON.stringify([placement.record_id,placement.peptide_id,placement.start,placement.end,role]);
    for(const e of cleavage.evidence||[]){
      if(e.role!=='spanning'&&e.is_candidate!==true)continue;
      const key=evidenceKey(e.site_id,e.record_id,e.role,e.peptide_id);
      if(!evidence.has(key))evidence.set(key,[]);
      evidence.get(key).push(e);
    }
    for(const t of cleavage.termini||[]){
      const key=terminalKey(t,t.role);
      if(!termini.has(key))termini.set(key,[]);
      termini.get(key).push(t);
    }
    const allCompatible=values=>!values.length?null:values.some(v=>v===false)?false:values.every(v=>v===true)?true:null;
    function compatible(placements,role){
      const values=[];
      for(const placement of placements){
        const matches=termini.get(terminalKey(placement,role))||[];
        if(!matches.length)values.push(null);
        for(const term of matches){
          if(!knownProfile||['reference_end','possible_initiator_met_removal'].includes(term.processing_status))values.push(null);
          else if(role==='N'&&term.pre_digest_n_label===true)values.push(false);
          else values.push(typeof term.digestion_compatible==='boolean'?term.digestion_compatible:null);
        }
      }
      return allCompatible(values);
    }
    function missedCleavages(id){
      const sequence=peptides.get(id)?.sequence;
      // Ambiguous residues can conceal or alter a cleavage context. U and O
      // are specified residues; X/B/Z/J are unresolved and cannot certify MC=0.
      if(!knownProfile||typeof sequence!=='string'||!sequence.length||!/^[ACDEFGHIKLMNPQRSTVWYOU]+$/.test(sequence))return null;
      let count=0;
      for(let i=0;i<sequence.length-1;i++)
        if((sequence[i]==='R'||(!profile.lysine_blocked&&sequence[i]==='K'))&&sequence[i+1]!=='P')count++;
      return count;
    }
    function otherCandidateBonds(pair,placements){
      const cacheKey=JSON.stringify([pair.site_id,pair.record_id,pair.spanning_peptide_id]);
      if(overlapCache.has(cacheKey))return overlapCache.get(cacheKey);
      const focal=sitesById.get(pair.site_id),positions=sitesByRecord.get(pair.record_id);
      if(!focal||focal.record_id!==pair.record_id||!Number.isInteger(focal.bond_left)||
        !positions||positions.some(k=>!Number.isInteger(k)||k<1)||!placements.length||
        placements.some(p=>!Number.isInteger(p.start)||!Number.isInteger(p.end)||p.start<1||p.end<=p.start||
          p.start>focal.bond_left||p.end<=focal.bond_left))return null;
      const found=new Set();
      for(const p of placements){
        let lo=0,hi=positions.length;
        while(lo<hi){const mid=(lo+hi)>>>1;if(positions[mid]<p.start)lo=mid+1;else hi=mid;}
        for(let i=lo;i<positions.length&&positions[i]<p.end;i++)if(positions[i]!==focal.bond_left)found.add(positions[i]);
      }
      overlapCache.set(cacheKey,found.size);return found.size;
    }
    return pair=>{
      const neo=evidence.get(evidenceKey(pair.site_id,pair.record_id,pair.role,pair.neo_peptide_id))||[];
      const span=evidence.get(evidenceKey(pair.site_id,pair.record_id,'spanning',pair.spanning_peptide_id))||[];
      return {neo_other_end_compatible:['N','C'].includes(pair.role)?compatible(neo,pair.role==='N'?'C':'N'):null,
        spanning_n_compatible:compatible(span,'N'),spanning_c_compatible:compatible(span,'C'),
        neo_missed_cleavages:missedCleavages(pair.neo_peptide_id),spanning_missed_cleavages:missedCleavages(pair.spanning_peptide_id),
        neo_placement_count:neo.length,spanning_placement_count:span.length,
        other_candidate_bonds:otherCandidateBonds(pair,span),
        profile_known:knownProfile};
    };
  }
  function selectionEligibility(pair,k){
    if(!pair.arithmetic_valid)return {eligible:false,reason:'invalid_arithmetic'};
    if(!pair.comparison)return {eligible:false,reason:'no_comparison'};
    const result=view(pair,k,1);
    const candidate=result.status==='quantified_in_both'||
      /^(neo|both)_quantified_only_in_(control|treatment)$/.test(result.status);
    return {eligible:candidate&&result.coverage_pass,
      reason:!candidate?result.status:result.coverage_pass?result.status:'below_selection_coverage'};
  }
  function rankPairs(pairs,options={}){
    const groups=new Set(pairs.map(p=>JSON.stringify([p.site_id,p.role,p.comparison])));
    if(groups.size>1)throw new Error('Representative selection requires one site, role and comparison.');
    const hasComparison=pairs.some(p=>Boolean(p.comparison));
    const inferred=pairs.length?(pairs[0].selection_k??(hasComparison?Math.max(1,Math.min(3,pairs[0].control.total,pairs[0].treatment.total)):0)):0;
    const k=options.selectionK??inferred;
    if(!Number.isInteger(k)||k<0)throw new Error('Selection coverage must be a non-negative integer.');
    if(hasComparison&&k===0)throw new Error('Selection coverage must be positive for a comparison.');
    if(options.selectionK===undefined&&pairs.some(p=>p.selection_k!==undefined&&p.selection_k!==k))
      throw new Error('Conflicting frozen selection coverage values.');
    const ranked=pairs.map(pair=>{
      const source=options.digestionByPair;
      const metadata=typeof source==='function'?source(pair):source instanceof Map?source.get(pair.key):source?.[pair.key];
      const declared=metadata??pair.digestion;
      const digestion=digestionPreference(declared);
      const overlap=declared?.other_candidate_bonds;
      const eligibility=selectionEligibility(pair,Math.max(1,k));
      return {...pair,selection:{rule_version:SELECTION_POLICY,selection_k:k,
        eligible:eligibility.eligible,eligibility_reason:eligibility.reason,
        preferred_digestion:digestion.preferred,digestion_reasons:digestion.reasons,
        other_candidate_bonds:Number.isInteger(overlap)&&overlap>=0?overlap:null,
        common_log2_intensity:finite(pair.common_log2_intensity)?pair.common_log2_intensity:null,
        fallback:false,rank:0}};
    });
    const fallback=!ranked.some(p=>p.selection.eligible);
    ranked.sort((a,b)=>{
      if(fallback)return coverageOrder(a,b)||identityOrder(a,b);
      const eligibleOrder=Number(b.selection.eligible)-Number(a.selection.eligible);
      if(eligibleOrder)return eligibleOrder;
      // Ineligible alternatives remain in legacy order for inspection.
      if(!a.selection.eligible)return coverageOrder(a,b)||identityOrder(a,b);
      const preference=Number(b.selection.preferred_digestion)-Number(a.selection.preferred_digestion);
      const intensityA=a.selection.common_log2_intensity,intensityB=b.selection.common_log2_intensity;
      const intensity=finite(intensityA)&&finite(intensityB)?intensityB-intensityA:
        Number(finite(intensityB))-Number(finite(intensityA));
      const overlapA=a.selection.other_candidate_bonds,overlapB=b.selection.other_candidate_bonds;
      const overlap=overlapA!==null&&overlapB!==null?overlapA-overlapB:
        Number(overlapA===null)-Number(overlapB===null);
      return preference||coverageOrder(a,b)||overlap||intensity||identityOrder(a,b);
    });
    ranked.forEach((pair,index)=>{pair.selection.fallback=fallback;pair.selection.rank=index+1;});
    return ranked;
  }
  function direction(value,tau){
    if(!finite(value))return null;
    if(tau===0)return Math.sign(value);
    return value>=tau?1:value<=-tau?-1:0;
  }
  const patterns={
    ratio_increased:'Ratio increased',ratio_decreased:'Ratio decreased',below_ratio_threshold:'Below ratio threshold',
    neo_up_spanning_down:'Neo up · spanning down',neo_down_spanning_up:'Neo down · spanning up',
    neo_up_spanning_below_threshold:'Neo up · spanning below threshold',neo_down_spanning_below_threshold:'Neo down · spanning below threshold',
    neo_below_threshold_spanning_up:'Neo below threshold · spanning up',neo_below_threshold_spanning_down:'Neo below threshold · spanning down',
    same_direction:'Both components change in the same direction',below_threshold:'Both components below threshold',
    only_in_one_condition:'Quantified in one condition',spanning_only:'Spanning quantified · neo not quantified',
    neo_only:'Neo quantified · spanning not quantified',insufficient_data:'Insufficient data',no_comparison:'No comparison'};
  function componentPattern(neo,span,tau){
    const a=direction(neo,tau),b=direction(span,tau);
    if(a===null||b===null)return 'insufficient_data';
    if(a===1&&b===-1)return 'neo_up_spanning_down';
    if(a===-1&&b===1)return 'neo_down_spanning_up';
    if(a===1&&b===0)return 'neo_up_spanning_below_threshold';
    if(a===-1&&b===0)return 'neo_down_spanning_below_threshold';
    if(a===0&&b===1)return 'neo_below_threshold_spanning_up';
    if(a===0&&b===-1)return 'neo_below_threshold_spanning_down';
    return a===0&&b===0?'below_threshold':'same_direction';
  }
  function view(pair,k=0,tau=1){
    if(!Number.isInteger(k)||k<0)throw new Error('Minimum coverage must be a non-negative integer.');
    if(!finite(tau)||tau<0)throw new Error('Effect threshold must be a non-negative finite number.');
    const c=pair.control,t=pair.treatment,threshold=Math.max(1,k);
    let status='insufficient_data',label='Insufficient data',pattern='insufficient_data',pass=false;
    const allMissing=(a,field)=>a.total>0&&a[field]===a.total;
    if(!pair.comparison){status='no_comparison';label='No comparison';pattern='no_comparison';}
    else if(pair.arithmetic_valid){
      if(finite(pair.log2_fc)){
        status='quantified_in_both';label='Quantified in both conditions';pass=c.co>=threshold&&t.co>=threshold;
        // The descriptive threshold applies to the within-sample ratio effect,
        // not separately to its two components. This does not change eligibility.
        const ratioDirection=direction(pair.log2_fc,tau);
        pattern=ratioDirection===1?'ratio_increased':ratioDirection===-1?'ratio_decreased':'below_ratio_threshold';
      }else if(t.co>0&&allMissing(c,'neoMissing')&&allMissing(c,'spanningMissing')){
        status='both_quantified_only_in_treatment';label='Both components quantified only in '+pair.treatment_name;pattern='only_in_one_condition';pass=t.co>=threshold;
      }else if(c.co>0&&allMissing(t,'neoMissing')&&allMissing(t,'spanningMissing')){
        status='both_quantified_only_in_control';label='Both components quantified only in '+pair.control_name;pattern='only_in_one_condition';pass=c.co>=threshold;
      }else if(t.co>0&&allMissing(c,'neoMissing')&&c.spanning>0){
        status='neo_quantified_only_in_treatment';label='Neo quantified only in '+pair.treatment_name;pattern='only_in_one_condition';pass=t.co>=threshold&&c.spanning>=threshold;
      }else if(c.co>0&&allMissing(t,'neoMissing')&&t.spanning>0){
        status='neo_quantified_only_in_control';label='Neo quantified only in '+pair.control_name;pattern='only_in_one_condition';pass=c.co>=threshold&&t.spanning>=threshold;
      }else if(c.neo+t.neo===0&&c.spanning+t.spanning>0){
        status='spanning_only';label=patterns.spanning_only;pattern='spanning_only';pass=Math.max(c.spanning,t.spanning)>=threshold;
      }else if(c.spanning+t.spanning===0&&c.neo+t.neo>0){
        status='neo_only';label=patterns.neo_only;pattern='neo_only';pass=Math.max(c.neo,t.neo)>=threshold;
      }
    }else{label='Inconsistent source values';}
    return {status,status_label:label,coverage_pass:k===0||pass,pattern,pattern_label:patterns[pattern]};
  }
  function mixedDirections(pairs,k=0,tau=1){
    const sites=new Set(pairs.map(p=>JSON.stringify([p.site_id,p.role,p.comparison])));
    if(sites.size>1)throw new Error('Mixed directions requires one site, role and comparison.');
    const eligible=pairs.filter(p=>view(p,k,tau).coverage_pass&&finite(p.log2_fc));
    return eligible.some(p=>direction(p.log2_fc,tau)===1)&&eligible.some(p=>direction(p.log2_fc,tau)===-1);
  }
  function directionSummary(pairs,k=0,tau=1){
    if(new Set(pairs.map(p=>JSON.stringify([p.role,p.comparison]))).size>1)throw new Error('Direction summary requires one role and comparison.');
    const unique=new Map();for(const pair of pairs){const key=measurementKey(pair);
      if(unique.has(key)){if(!close(unique.get(key).log2_fc,pair.log2_fc))throw new Error('Conflicting effects for one measurement pair.');}
      else unique.set(key,pair);
    }
    const eligible=[...unique.values()].filter(p=>p.arithmetic_valid&&view(p,k,tau).coverage_pass&&finite(p.log2_fc));
    const increased=eligible.filter(p=>direction(p.log2_fc,tau)===1).length,decreased=eligible.filter(p=>direction(p.log2_fc,tau)===-1).length;
    const directional=increased+decreased,below_threshold=eligible.length-directional;
    let label='Below threshold';
    if(!eligible.length)label='No computable pairs';
    else if(increased&&!decreased)label='Increased';else if(decreased&&!increased)label='Decreased';
    else if(directional){const fraction=Math.max(increased,decreased)/directional;label=fraction>=.8?(increased>decreased?'Predominantly increased':'Predominantly decreased'):'Mixed directions';}
    return {label,increased,decreased,below_threshold,computable:eligible.length,total:unique.size,directional};
  }
  function groupOverviewRows(rows){
    const groups=new Map();
    for(const row of rows){const p=row.selected,key=p?'measurement:'+measurementKey(p):'mapping:'+row.key;
      if(!groups.has(key))groups.set(key,{...row,key,mappings:[]});
      const group=groups.get(key),first=group.selected;
      if(first&&p){
        const fields=['log2_fc','neo_log2_fc','spanning_log2_fc','common_log2_intensity'];
        if(fields.some(f=>!close(first[f],p[f]))||['control','treatment'].some(side=>['co','total','neo','spanning'].some(f=>first[side][f]!==p[side][f]))||first.arithmetic_valid!==p.arithmetic_valid)
          throw new Error('Alternative mappings have inconsistent measurement-pair quantities.');
      }
      group.mappings.push(row);
    }
    for(const group of groups.values()){
      const pairs=new Map();for(const row of group.mappings)for(const pair of row.pairs){const key=measurementKey(pair);if(!pairs.has(key))pairs.set(key,pair);}
      group.pairs=[...pairs.values()];group.computable=group.pairs.filter(p=>finite(p.log2_fc)).length;
    }
    return [...groups.values()];
  }
  // Descriptive residue composition. No enrichment test or p-value is calculated.
  // Callers pass frozen representatives, never a minimum-p alternative.
  function motifSummary(rows,{role,family='',k=0,tau=1}={}){
    if(!['N','C'].includes(role))throw new Error('Motif summary requires one terminal role.');
    if(typeof family!=='string')throw new Error('Motif summary requires one test-family label.');
    if(!Number.isInteger(k)||k<0||!finite(tau)||tau<0)throw new Error('Invalid motif coverage or ratio threshold.');
    const residues='ACDEFGHIKLMNPQRSTVWY'.split(''),positions=['P4','P3','P2','P1','P1prime','P2prime','P3prime','P4prime'];
    const unique=new Map();
    for(const row of rows||[]){
      const pair=row.pair,site=row.site,families=[...new Set(row.families||[])];
      if(!pair||!site||pair.role!==role||!finite(pair.log2_fc)||!pair.arithmetic_valid||!view(pair,k,tau).coverage_pass)continue;
      if(family? !families.includes(family):families.length)continue;
      if(site.record_id!==pair.record_id||site.site_id!==pair.site_id||!Number.isInteger(site.bond_left))
        throw new Error('Motif representative does not match its reference bond.');
      if(['comparison','neo_peptide_id','spanning_peptide_id'].some(field=>typeof pair[field]!=='string'||!pair[field]))
        throw new Error('Motif representatives require comparison and peptide-form identifiers.');
      const key=measurementKey(pair),mappingKey=JSON.stringify([site.record_id,site.bond_left]);
      if(!unique.has(key))unique.set(key,{pair,mappings:new Map()});
      const group=unique.get(key);
      if(!close(group.pair.log2_fc,pair.log2_fc)||['control','treatment'].some(side=>
        ['co','total'].some(field=>group.pair[side][field]!==pair[side][field])))
        throw new Error('Conflicting motif quantities for one measurement pair.');
      const existing=group.mappings.get(mappingKey);
      if(existing&&positions.some(p=>existing[p]!==site[p]))throw new Error('Conflicting motif contexts for one protein-record/bond.');
      group.mappings.set(mappingKey,site);
    }
    // Alternative mappings are not additional measurements. Retain a residue
    // only when all mappings agree at that position; do not pick a protein or
    // distribute fractional counts across incompatible reference contexts.
    const measured=[...unique.values()].map(group=>({pair:group.pair,site:Object.fromEntries(positions.map(position=>{
      const values=[...group.mappings.values()].map(site=>site[position]);
      return [position,residues.includes(values[0])&&values.every(value=>value===values[0])?values[0]:null];
    }))}));
    const up=measured.filter(r=>direction(r.pair.log2_fc,tau)===1),
      down=measured.filter(r=>direction(r.pair.log2_fc,tau)===-1);
    const count=group=>positions.map(position=>{
      const valid=group.map(r=>r.site[position]).filter(r=>residues.includes(r));
      return {denominator:valid.length,counts:Object.fromEntries(residues.map(r=>[r,valid.filter(x=>x===r).length]))};
    });
    const make=(name,group)=>{
      const counts=count(group),cells=[];
      for(const residue of residues)for(const [i,position] of positions.entries()){
        const numerator=counts[i].counts[residue],denominator=counts[i].denominator;
        cells.push({group:name,position,residue,numerator,denominator,percent:denominator?100*numerator/denominator:null});
      }
      return {name,measurement_pairs:group.length,cells};
    };
    return {role,family,coverage_k:k,threshold:tau,measurement_pairs:measured.length,residues,positions,
      groups:[make('Ratio increased',up),make('Ratio decreased',down)]};
  }
  return Object.freeze({SELECTION_POLICY,pairKey,measurementKey,selectedSamples,defaultK,summarizePairs,preparedPairs,preparedOrder,rankPairs,legacyRankPairs,
    digestionPreference,buildDigestionLookup,selectionEligibility,view,direction,componentPattern,mixedDirections,directionSummary,groupOverviewRows,motifSummary,patterns});
});
