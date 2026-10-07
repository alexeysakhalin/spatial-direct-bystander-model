"""Reconstruct spatial observables from the central saved simulations."""
from pathlib import Path
import argparse,hashlib,json,platform
import numpy,pandas,scipy
import geometry,gradients
from ratio_metrics import metric_values,units
from summaries import descriptive,difference
from projection import project
HERE=Path(__file__).resolve().parent

def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):
    if p.exists():raise ValueError('Output already exists: '+p.name)
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-root',type=Path,required=True)
    p.add_argument('--inputs-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    assert not a.output.exists(),'Choose a new output directory'
    for n,h in load(HERE/'MANIFEST_SHA256.json').items():assert sha(HERE/n)==h,n
    technical=load(HERE/'TECHNICAL_ACCEPTANCE.json');cases=load(HERE/'CASES.json');expected=load(HERE/'EXPECTED_SPATIAL_ROWS.json');counts=load(HERE/'EXPECTED_COUNTS.json');raw=load(HERE/'RAW_FILES.json')
    assert len(cases)==18 and {r['id']for r in cases}==set(expected)==set(raw)==set(counts)
    assert {(r['ratio'],r['arm'],r['seed'])for r in cases}=={(r,arm,s)for r in[2,5]for arm in['2d','3d','local_entry']for s in range(3)}
    tests=geometry.selftest()+gradients.selftests();geometry.REGISTRY=HERE/'CELL_STATES.json'
    a.output.mkdir(parents=True);values={};case_metrics=[];files=0;rawbytes=0;hourly=0
    for c in cases:
        cid=c['id'];out=a.raw_root.resolve()/cid/'output';cfg=a.inputs_root.resolve()/cid
        assert len(list(out.glob('output[0-9]*.xml')))==97
        for n,h in c['portable_input_sha256'].items():assert sha(cfg/n)==h
        for n,r in raw[cid].items():
            f=out/n;assert Path(n).name==n and f.is_file()and not f.is_symlink()and f.stat().st_size==r['bytes']and sha(f)==r['sha256'],cid+'/'+n
            files+=1;rawbytes+=r['bytes']
        rec=dict(c,output=str(out),folder=str(cfg),input_sha256=c['portable_input_sha256'])
        accepted=technical[cid];assert accepted['last']['live_targets']==counts[cid][-1]['live_targets']
        ar=geometry.inspect(rec,accepted);gr=gradients.inspect(rec,accepted,HERE/'CELL_STATES.json')
        # The accepted-record fields refer to the bound prior technical validation.
        # This invocation recomputes geometry and saved-field derivatives.
        actual=project(c,ar,gr);assert actual==expected[cid],cid
        used=dict(ar['raw_cell_input_sha256']);used.update(gr['raw_sha256'])
        assert {Path(n).name for n in used}==set(raw[cid])
        for n,h in used.items():assert h==raw[cid][Path(n).name]['sha256'],cid+'/'+Path(n).name
        for h in range(97):
            assert ar['rows'][h]['time_h']==gr['rows'][h]['time_h']==h
            assert ar['rows'][h]['living_targets']==gr['rows'][h]['counts']['live_targets']==counts[cid][h]['live_targets']
            assert ar['rows'][h]['living_effectors']==gr['rows'][h]['counts']['live_effectors']==counts[cid][h]['live_effectors']
            assert gr['rows'][h]['counts']['retained_dead']==counts[cid][h]['retained_dead'];hourly+=1
        save(a.output/'cases'/f'{cid}.json',actual)
        m=metric_values(ar,gr);assert len(m)==59;values[c['ratio'],c['arm'],c['seed']]=m
        case_metrics.append({'case':cid,'ratio':c['ratio'],'arm':c['arm'],'seed':c['seed'],'metrics':m})
    summaries=[];pairs=[];pair_summaries=[];keys=list(next(iter(values.values())))
    for ratio in[2,5]:
        for arm in['2d','3d','local_entry']:
            for key in keys:summaries.append({'ratio':ratio,'arm':arm,'metric':key,'unit':units(key),'seeds':[0,1,2],**descriptive([values[ratio,arm,s][key]for s in range(3)])})
    for arm in['2d','3d','local_entry']:
        for key in keys:
            ds=[]
            for seed in range(3):
                d=difference(values[2,arm,seed][key],values[5,arm,seed][key]);pairs.append({'arm':arm,'seed':seed,'metric':key,'unit':units(key),**d});ds.append(d['TE5_minus_TE2'])
            pair_summaries.append({'arm':arm,'metric':key,'unit':units(key),'same_seed_differences':ds,**descriptive(ds)})
    assert len(summaries)==354 and len(pairs)==531 and len(pair_summaries)==177 and hourly==1746 and files==5274
    for name,data in [('N3_SPATIAL_METRIC_SUMMARIES',summaries),('PAIRED_RATIO_SPATIAL_DIFFERENCES',pairs),('N3_RATIO_SPATIAL_DIFFERENCES',pair_summaries)]:
        assert data==load(HERE/('EXPECTED_'+name+'.json')),name;save(a.output/(name+'.json'),data)
    save(a.output/'CASE_METRICS.json',case_metrics)
    save(a.output/'ENVIRONMENT.json',{'python':platform.python_version(),'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__})
    save(a.output/'REPLAY_VERIFICATION.json',{'passed':True,'central_cases':18,'hourly_access_states':1746,'hourly_CXCL_field_states':1746,'initial_target_hulls':18,'case_metric_values':1062,'N3_scenario_metric_summaries':354,'paired_ratio_metric_differences':531,'paired_ratio_N3_summaries':177,'consumed_raw_files_hash_verified':files,'consumed_raw_bytes_hash_verified':rawbytes,'portable_input_files_verified':54,'geometry_and_gradient_selftests':tests,'prior_technical_acceptance_bound':True,'new_simulations':0,'numerical_convergence_established':False,'experimental_validation_established':False})
    save(a.output/'RESULT_MANIFEST_SHA256.json',{str(f.relative_to(a.output)):sha(f)for f in sorted(a.output.rglob('*'))if f.is_file()})
    print(json.dumps({'passed':True,'central_cases':18,'hourly_states_per_observable':hourly,'raw_files_verified':files}))
if __name__=='__main__':main()
