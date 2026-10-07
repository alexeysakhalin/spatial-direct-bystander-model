"""Reconstruct accepted controls and numerical checks from saved outputs."""
from pathlib import Path
import argparse,json,math,platform,sys,xml.etree.ElementTree as ET
import numpy,pandas,scipy
from read_counts import digest,read_cells
from saved_fields import initial_audit,audit_frame
from scalar_metrics import enrich_counts,metric_values,selftest
from medium import medium_audit
HERE=Path(__file__).resolve().parent

def load(p):return json.loads(p.read_text())
def save(p,d):
 assert not p.exists(),p.name
 p.write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')
def close(a,b):
 if isinstance(b,dict):
  assert set(a)==set(b)
  for k in b:close(a[k],b[k])
 elif isinstance(b,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):close(x,y)
 elif isinstance(b,float):assert a is not None and math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)
 else:assert a==b,(a,b)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--raw-root',required=True,type=Path);ap.add_argument('--inputs-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();assert not args.output.exists(),'Choose a new output directory'
 for n,h in load(HERE/'MANIFEST_SHA256.json').items():assert digest(HERE/n)==h,n
 records=load(HERE/'CASES.json');payload=load(HERE/'RAW_FILES.json');assert len(records)==5 and set(payload)=={r['id']for r in records}
 registry=HERE/'CELL_STATES.json';reg={x['id']:x for x in load(registry)['states']};base={(x['ratio'],x['arm'],x['seed']):x for x in load(HERE/'CENTRAL_METRICS.json')};assert len(base)==18
 tests=selftest();assert len(tests)==12;args.output.mkdir(parents=True);(args.output/'cases').mkdir();allmetrics=[];differences=[];rawfiles=rawbytes=full=partial=0;maxresidual=0.;summaries=[]
 for r in records:
  cid=r['id'];out=args.raw_root/cid/'output';cfg=args.inputs_root/cid;ref=load(HERE/'references'/f'{cid}.json');assert r['technical_exit_code']==ref['exit_code']==0
  for n,h in r['portable_input_sha256'].items():assert digest(cfg/n)==h,cid+'/'+n
  assert len(list(out.glob('output[0-9]*.xml')))==97
  for n,x in payload[cid].items():
   assert Path(n).name==n;p=out/n;assert p.is_file()and not p.is_symlink()and p.stat().st_size==x['bytes']and digest(p)==x['sha256'],cid+'/'+n
   rawfiles+=1;rawbytes+=x['bytes']
  rec=dict(r,folder=str(cfg));initial=initial_audit(out,rec,registry)
  xmls=[out/'initial.xml']+[out/f'output{h:08d}.xml'for h in range(97)]+[out/'final.xml'];rows=[audit_frame(p,registry,r['no_effectors_control'])for p in xmls]
  assert len(rows)==99 and rows[0]['time_h']==0;dt=float(ET.parse(cfg/'settings.xml').findtext('overall/dt_diffusion'));assert abs(rows[-1]['time_h']-96-dt/60)<1e-7
  counts=[read_cells(p,reg)for p in xmls[1:-1]];assert counts==ref['counts'];hourly=enrich_counts(rows[1:-1],counts);assert all(x['retained_dead']==y['retained_dead']for x,y in zip(hourly,counts))
  if ref['fields']is not None:close(rows,ref['fields']);full+=1
  else:partial+=1
  m=medium_audit(out,rec);assert m['passed']and m['numbered_frames']==97 and all(x['ledger_rows']==98 for x in m['species'].values())
  if ref['medium']['validator_version']==3:close(m,ref['medium'])
  else:
   original=ref['medium']['field_time_series'];assert len(original)==98
   for actual,prior in zip(m['field_time_series'],original):
    assert actual['snapshot']==prior['snapshot'];close(actual['time_min']/60,prior['time_h'])
    for sp in['ifng','tnf']:close(actual[sp],prior[sp])
  maxresidual=max(maxresidual,*(x['max_independent_relative_residual']for x in m['species'].values()))
  vals={k:v[0]for k,v in metric_values(hourly).items()};assert len(vals)==150;central=base[r['reference_ratio'],r['arm'],r['seed']];assert set(vals)==set(central['values'])
  allmetrics.append({'case':cid,'arm':r['arm'],'seed':r['seed'],'target_to_effector_ratio':r['target_to_effector_ratio'],'values':vals})
  for key,v in vals.items():
   refv=central['values'][key];differences.append({'case':cid,'reference_case':central['case'],'reference_ratio':r['reference_ratio'],'no_effectors_control':r['no_effectors_control'],'metric':key,'value':v,'reference_value':refv,'difference':None if v is None or refv is None else v-refv})
  save(args.output/'cases'/f'{cid}.json',{'case':cid,'arm':r['arm'],'seed':r['seed'],'target_to_effector_ratio':r['target_to_effector_ratio'],'initial':initial,'initial_field_state':rows[0],'hourly_rows':hourly,'terminal_field_state':rows[-1],'medium_balance':m,'native_full_field_reference':r['native_full_field_reference']})
  summaries.append({'case':cid,'native_exit_code':0,'hourly_count_states':97,'all_field_states':99,'balance_rows':196,'native_full_field_reference':r['native_full_field_reference'],'max_independent_relative_balance_residual':max(x['max_independent_relative_residual']for x in m['species'].values())})
  print(json.dumps({'case':cid,'hourly_states':97,'field_states':99,'balance_rows':196,'native_full_field_reference':r['native_full_field_reference']}),flush=True)
 assert full==5 and partial==0 and len(differences)==750
 save(args.output/'CASE_METRICS.json',allmetrics);save(args.output/'PAIRED_CENTRAL_DIFFERENCES.json',differences);save(args.output/'TECHNICAL_REPLAY_SUMMARY.json',summaries)
 save(args.output/'REPLAY_VERIFICATION.json',{'passed':True,'accepted_noncentral_cases':5,'hourly_population_states':485,'all_field_states':495,'fields_per_state':6,'balance_rows':980,'case_metric_values':750,'paired_central_metric_differences':750,'native_full_field_reference_cases':full,'native_partial_field_reference_cases':partial,'consumed_raw_files_hash_verified':rawfiles,'consumed_raw_bytes_hash_verified':rawbytes,'portable_input_files_hash_verified':15,'max_independent_relative_balance_residual':maxresidual,'scientific_boundary_tests':tests,'new_simulations':0,'numerical_convergence_accepted':False,'biological_calibration_established':False,'scope':'Replay of accepted saved checks. No new materiality or convergence thresholds are applied. All five extension cases have complete native references.'})
 save(args.output/'ENVIRONMENT.json',{'python':sys.version,'implementation':platform.python_implementation(),'machine':platform.machine(),'byteorder':sys.byteorder,'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__})
 save(args.output/'RESULT_MANIFEST_SHA256.json',{str(p.relative_to(args.output)):digest(p)for p in sorted(args.output.rglob('*'))if p.is_file()})
if __name__=='__main__':main()
