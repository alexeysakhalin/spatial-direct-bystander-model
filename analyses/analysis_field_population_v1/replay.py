"""Reproduce population and field summaries for the central simulations."""
from pathlib import Path
import argparse,json,math,platform,sys,xml.etree.ElementTree as ET
import numpy,pandas,scipy
from read_counts import digest,read_cells
from saved_fields import initial_audit,audit_frame
from scalar_metrics import enrich_counts,metric_values,selftest
from ensemble_stats import stats
HERE=Path(__file__).resolve().parent

def load(p):return json.loads(p.read_text())
def save(p,d):
    assert not p.exists(),p.name
    p.write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')
def compare(a,b):
    if isinstance(b,dict):
        assert set(a)==set(b)
        for k in b:compare(a[k],b[k])
    elif isinstance(b,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(b,float):assert a is not None and math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)
    else:assert a==b,(a,b)
def unit(key):
    if key.startswith('count/'):return 'cells'
    return 'relative field units * hours' if key.endswith('/AUC_0_96h') else 'relative field units'
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--raw-root',type=Path,required=True);ap.add_argument('--inputs-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert not args.output.exists(),'Choose a new output directory'
    for n,h in load(HERE/'MANIFEST_SHA256.json').items():assert digest(HERE/n)==h,n
    cases=load(HERE/'CASES.json');payload=load(HERE/'RAW_FILES.json');expected=load(HERE/'EXPECTED_HOURLY_ROWS.json');counts=load(HERE/'EXPECTED_COUNTS.json');accepted=load(HERE/'TECHNICAL_ACCEPTANCE.json')
    assert len(cases)==18 and set(payload)==set(expected)==set(counts)==set(accepted)=={r['id']for r in cases}
    assert {(r['ratio'],r['arm'],r['seed'])for r in cases}=={(ratio,arm,seed)for ratio in[2,5]for arm in['2d','3d','local_entry']for seed in range(3)}
    registry=HERE/'CELL_STATES.json';reg={r['id']:r for r in load(registry)['states']};tests=selftest();assert len(tests)==12
    args.output.mkdir(parents=True);(args.output/'cases').mkdir();metrics=[];rawfiles=rawbytes=0
    for rec in cases:
        cid=rec['id'];root=args.raw_root/cid/'output';cfg=args.inputs_root/cid;acc=accepted[cid]
        assert acc['exit_code']==0 and acc['frames']==97 and acc['end_h']==96 and acc['medium_audit_passed']
        assert len(list(root.glob('output[0-9]*.xml')))==97
        for n,h in rec['portable_input_sha256'].items():assert digest(cfg/n)==h
        for n,x in payload[cid].items():
            assert Path(n).name==n;p=root/n;assert p.is_file()and not p.is_symlink()and p.stat().st_size==x['bytes']and digest(p)==x['sha256'],cid+'/'+n
            rawfiles+=1;rawbytes+=x['bytes']
        initial=initial_audit(root,dict(rec,folder=str(cfg)),registry)
        xmls=[root/'initial.xml']+[root/f'output{h:08d}.xml'for h in range(97)]+[root/'final.xml']
        rows=[audit_frame(p,registry)for p in xmls];assert len(rows)==99 and rows[0]['time_h']==0
        dt=float(ET.parse(cfg/'settings.xml').findtext('overall/dt_diffusion'));assert abs(rows[-1]['time_h']-(96+dt/60))<1e-7
        actual_counts=[read_cells(p,reg)for p in xmls[1:-1]];assert actual_counts==counts[cid]
        hourly=enrich_counts(rows[1:-1],actual_counts);compare(hourly,expected[cid])
        assert all(r['retained_dead']==c['retained_dead']for r,c in zip(hourly,actual_counts))
        values=metric_values(hourly);assert len(values)==150
        metrics.append({'case':cid,'ratio':rec['ratio'],'arm':rec['arm'],'seed':rec['seed'],'values':{k:v[0]for k,v in values.items()}})
        save(args.output/'cases'/f'{cid}.json',{'case':cid,'ratio':rec['ratio'],'arm':rec['arm'],'seed':rec['seed'],'initial':initial,'initial_field_state':rows[0],'hourly_rows':hourly,'terminal_field_state':rows[-1]})
        print(json.dumps({'case':cid,'hourly_states':97,'all_field_states':99,'metrics':150}),flush=True)
    groups=[];differences=[];difference_groups=[]
    by={(r['ratio'],r['arm'],r['seed']):r for r in metrics};keys=sorted(metrics[0]['values']);assert all(set(m['values'])==set(keys)for m in metrics)
    for ratio in[2,5]:
        for arm in['2d','3d','local_entry']:
            for k in keys:
                vv=[by[ratio,arm,s]['values'][k]for s in range(3)];groups.append({'ratio':ratio,'arm':arm,'metric':k,'units':unit(k),'seed_values':vv,**stats(vv)})
    contrasts=[]
    for ratio in[2,5]:
        for left,right in[('3d','2d'),('local_entry','2d'),('local_entry','3d')]:contrasts.append((f'TE{ratio}_{left}_minus_{right}',ratio,left,ratio,right))
    for arm in['2d','3d','local_entry']:contrasts.append((f'{arm}_TE5_minus_TE2',5,arm,2,arm))
    for name,lr,la,rr,ra in contrasts:
        for k in keys:
            deltas=[]
            for seed in range(3):
                left=by[lr,la,seed]['values'][k];right=by[rr,ra,seed]['values'][k];delta=left-right if left is not None and right is not None else None;deltas.append(delta)
                differences.append({'contrast':name,'left_ratio':lr,'left_arm':la,'right_ratio':rr,'right_arm':ra,'seed':seed,'metric':k,'units':unit(k),'left':left,'right':right,'difference':delta})
            difference_groups.append({'contrast':name,'metric':k,'units':unit(k),'seed_differences':deltas,**stats(deltas)})
    assert len(groups)==900 and len(differences)==4050 and len(difference_groups)==1350
    save(args.output/'CASE_METRICS.json',metrics);save(args.output/'N3_METRIC_SUMMARIES.json',groups);save(args.output/'PAIRED_DIFFERENCES.json',differences);save(args.output/'N3_DIFFERENCE_SUMMARIES.json',difference_groups)
    save(args.output/'REPLAY_VERIFICATION.json',{'passed':True,'central_cases':18,'hourly_population_states':1746,'all_field_states':1782,'fields_per_state':6,'metrics_per_case':150,'case_metric_values':2700,'N3_metric_summaries':900,'paired_metric_differences':4050,'N3_difference_summaries':1350,'consumed_raw_files_hash_verified':rawfiles,'consumed_raw_bytes_hash_verified':rawbytes,'portable_input_files_hash_verified':54,'qualified_hourly_rows_reproduced':1746,'scientific_boundary_tests':tests,'new_simulations':0,'scope':'Population and field summaries of central saved outputs. Three-seed ranges describe these realizations; they are not confidence intervals. This replay does not establish numerical convergence or biological calibration.'})
    save(args.output/'ENVIRONMENT.json',{'python':sys.version,'implementation':platform.python_implementation(),'machine':platform.machine(),'byteorder':sys.byteorder,'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__})
    save(args.output/'RESULT_MANIFEST_SHA256.json',{str(p.relative_to(args.output)):digest(p)for p in sorted(args.output.rglob('*'))if p.is_file()})
if __name__=='__main__':main()
