"""Reconstruct spatial observables of accepted model checks."""
from pathlib import Path
import argparse,hashlib,json,math,platform
import numpy,pandas,scipy
import geometry,gradients
from ratio_metrics import metric_values
from projection import project
HERE=Path(__file__).resolve().parent

def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):
 assert not p.exists(),p.name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--raw-root',type=Path,required=True);ap.add_argument('--inputs-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists(),'Choose a new output directory'
 for n,h in load(HERE/'MANIFEST_SHA256.json').items():assert sha(HERE/n)==h,n
 cases=load(HERE/'CASES.json');raw=load(HERE/'RAW_FILES.json');assert len(cases)==2 and set(raw)=={c['id']for c in cases}
 central={(c['ratio'],c['arm'],c['seed']):c for c in load(HERE/'CENTRAL_METRICS.json')};assert len(central)==18
 tests=geometry.selftest()+gradients.selftests();geometry.REGISTRY=HERE/'CELL_STATES.json';a.output.mkdir(parents=True);cms=[];diffs=[];files=rawbytes=hourly=reference_cases=0
 for c in cases:
  cid=c['id'];out=a.raw_root.resolve()/cid/'output';cfg=a.inputs_root.resolve()/cid;ref=load(HERE/'references'/f'{cid}.json');counts=ref['counts'];accepted=ref['acceptance'];assert accepted['exit_code']==0 and accepted['frames']==97 and accepted['end_h']==96 and accepted['medium_audit_passed']
  for n,h in c['portable_input_sha256'].items():assert sha(cfg/n)==h,cid+'/'+n
  assert len(list(out.glob('output[0-9]*.xml')))==97
  for n,r in raw[cid].items():
   f=out/n;assert Path(n).name==n and f.is_file()and not f.is_symlink()and f.stat().st_size==r['bytes']and sha(f)==r['sha256'],cid+'/'+n;files+=1;rawbytes+=r['bytes']
  rec=dict(c,output=str(out),folder=str(cfg),input_sha256=c['portable_input_sha256']);ar=geometry.inspect(rec,accepted);gr=gradients.inspect(rec,accepted,HERE/'CELL_STATES.json');actual=project(c,ar,gr)
  if ref['native_spatial_rows']is not None:assert actual==ref['native_spatial_rows'],cid;reference_cases+=1
  used=dict(ar['raw_cell_input_sha256']);used.update(gr['raw_sha256']);assert {Path(n).name for n in used}==set(raw[cid])
  for n,h in used.items():assert h==raw[cid][Path(n).name]['sha256'],cid+'/'+Path(n).name
  for h in range(97):
   aa=ar['rows'][h];gg=gr['rows'][h];cc=counts[h];assert aa['time_h']==gg['time_h']==cc['time_h']==h
   assert aa['living_targets']==gg['counts']['live_targets']==cc['live_targets']and aa['living_effectors']==gg['counts']['live_effectors']==cc['live_effectors']and gg['counts']['retained_dead']==cc['retained_dead'];hourly+=1
  m=metric_values(ar,gr);assert len(m)==59
  if ref['native_spatial_metrics']is not None:
   for k,v in ref['native_spatial_metrics'].items():assert m[k]==v,(cid,k,m[k],v)
  cms.append({'case':cid,'ratio':c['ratio'],'reference_ratio':c['reference_ratio'],'arm':c['arm'],'seed':c['seed'],'metrics':m,'native_spatial_reference':c['native_spatial_reference']});base=central[c['reference_ratio'],c['arm'],c['seed']]
  for k,v in m.items():
   bv=base['metrics'][k];diffs.append({'case':cid,'reference_case':base['case'],'metric':k,'value':v,'reference_value':bv,'difference':None if v is None or bv is None else v-bv})
  save(a.output/'cases'/f'{cid}.json',dict(actual,native_spatial_reference=c['native_spatial_reference']));print(json.dumps({'case':cid,'hourly_states':97,'native_spatial_reference':c['native_spatial_reference']}),flush=True)
 assert files==586 and hourly==194 and reference_cases==2 and len(diffs)==118
 save(a.output/'CASE_METRICS.json',cms);save(a.output/'PAIRED_CENTRAL_DIFFERENCES.json',diffs)
 save(a.output/'REPLAY_VERIFICATION.json',{'passed':True,'accepted_noncentral_cases':2,'hourly_access_states':194,'hourly_CXCL_field_states':194,'case_metric_values':118,'paired_central_metric_differences':118,'native_spatial_reference_cases':2,'newly_derived_medium_comparison_cases':0,'consumed_raw_files_hash_verified':files,'consumed_raw_bytes_hash_verified':rawbytes,'portable_input_files_verified':6,'geometry_and_gradient_selftests':tests,'new_simulations':0,'numerical_convergence_established':False,'experimental_validation_established':False})
 save(a.output/'ENVIRONMENT.json',{'python':platform.python_version(),'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__});save(a.output/'RESULT_MANIFEST_SHA256.json',{str(f.relative_to(a.output)):sha(f)for f in sorted(a.output.rglob('*'))if f.is_file()})
 print(json.dumps({'passed':True,'cases':2,'hourly_states':hourly,'raw_files_verified':files}))
if __name__=='__main__':main()
