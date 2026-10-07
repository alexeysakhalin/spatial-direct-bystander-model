"""Reproduce central population trajectories and medium balances from saved data."""
from pathlib import Path
import argparse,csv,json,math,platform,sys,xml.etree.ElementTree as ET
from read_saved_data import digest,read_cells,read_field_v4

HERE=Path(__file__).resolve().parent

def read(path):
    return json.loads(path.read_text())

def write(path,value):
    if path.exists():raise ValueError('Output already exists: '+path.name)
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def integral(values):
    return math.fsum((a+b)/2 for a,b in zip(values,values[1:]))

def stats(values):
    return {'mean':math.fsum(values)/len(values),'min':min(values),'max':max(values)}

def close_nested(actual,expected):
    if isinstance(expected,dict):
        assert set(actual)==set(expected)
        for k in expected:close_nested(actual[k],expected[k])
    elif isinstance(expected,list):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):close_nested(a,b)
    elif isinstance(expected,(int,float)):
        assert math.isfinite(actual) and abs(actual-expected)<=1e-10*(1+abs(expected)),(actual,expected)
    else:assert actual==expected,(actual,expected)

def verify_balance(rec,raw,settings):
    cfg=ET.parse(settings);up=cfg.find('user_parameters')
    reservoir=float(up.findtext('ikm_medium_reservoir_volume_um3'))
    flow=float(up.findtext('ikm_medium_external_flow_um3_per_min'));dt=float(cfg.findtext('overall/dt_diffusion'))
    ti,vi,initial=read_field_v4(raw/'initial.xml')
    assert abs(ti)<1e-7 and abs(vi-rec['domain_volume_um3'])<1e-8*vi
    ledger=list(csv.DictReader((raw/'medium_balance.csv').open()))
    by={k:[{n:float(v)for n,v in r.items()if n!='substrate'}for r in ledger if r['substrate']==k]for k in ['ifng','tnf']}
    assert len(ledger)==196 and all(len(v)==98 for v in by.values())
    for key,rows in by.items():
        initial[key]+=reservoir*float(up.findtext('ikm_medium_'+key+'_initial'))
        assert all(math.isfinite(x)for r in rows for x in r.values())
        assert all(b['time_min']>a['time_min']for a,b in zip(rows,rows[1:]))
    maximum=0.
    for j,p in enumerate([raw/f'output{i:08d}.xml'for i in range(97)]+[raw/'final.xml']):
        t,volume,masses=read_field_v4(p)
        assert abs(t-(60*j if j<97 else 5760+dt))<1e-6
        assert abs(volume-rec['domain_volume_um3'])<1e-8*volume
        for key,mass in masses.items():
            r=by[key][j];assert abs(r['time_min']-t)<1e-6
            assert abs(r['domain_amount']-mass)<1e-8*(1+abs(mass))
            assert abs(r['reservoir_amount']-reservoir*r['reservoir_concentration'])<1e-8*(1+abs(r['reservoir_amount']))
            assert r['reservoir_decay_loss']>=0 and r['reservoir_concentration']>=0
            if flow==0:assert r['external_net_export']==0
            residual=mass+r['reservoir_amount']-initial[key]-r['native_net_change']+r['external_net_export']+r['reservoir_decay_loss']
            scale=1+abs(initial[key])+abs(r['native_net_change'])+abs(r['external_net_export'])+abs(r['reservoir_decay_loss'])
            assert abs(residual)<1e-8*scale and abs(r['balance_residual']-residual)<1e-8*scale
            maximum=max(maximum,abs(residual)/scale)
    return {'case':rec['id'],'field_states':99,'ledger_rows':196,'max_relative_balance_residual':maximum,'balanced_substrates':['ifng','tnf']}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--raw-root',required=True,type=Path,help='Directory containing one subdirectory per case')
    ap.add_argument('--inputs-root',required=True,type=Path,help='Directory containing the portable case configurations')
    ap.add_argument('--output',required=True,type=Path,help='New directory for derived results')
    args=ap.parse_args();rawroot=args.raw_root.resolve();inputs=args.inputs_root.resolve()
    if args.output.exists():raise ValueError('Choose a new output directory')
    for name,h in read(HERE/'MANIFEST_SHA256.json').items():assert digest(HERE/name)==h,name
    cases=read(HERE/'CASES.json');assert len(cases)==18
    assert {(r['ratio'],r['arm'],r['seed'])for r in cases}=={(r,a,s)for r in [2,5]for a in ['2d','3d','local_entry']for s in range(3)}
    registry={r['id']:r for r in read(HERE/'CELL_STATES.json')['states']};byname={r['name']:r for r in registry.values()}
    expected=read(HERE/'EXPECTED_COUNTS.json');payload=read(HERE/'RAW_FILES.json')
    assert set(expected)==set(payload)=={r['id']for r in cases}
    args.output.mkdir(parents=True)
    hourly=[];metrics=[];balances=[];rawfiles=0;rawbytes=0;trajectories={}
    for rec in cases:
        cid=rec['id'];out=rawroot/cid/'output';cfg=inputs/cid
        assert len(list(out.glob('output[0-9]*.xml')))==97
        for n,h in rec['portable_input_sha256'].items():assert digest(cfg/n)==h,cid+'/'+n
        assert set(payload[cid])=={p.name for p in out.iterdir()if p.is_file()}
        for n,r in payload[cid].items():
            assert Path(n).name==n;p=out/n
            assert p.is_file()and not p.is_symlink()and p.stat().st_size==r['bytes']and digest(p)==r['sha256'],cid+'/'+n
            rawfiles+=1;rawbytes+=r['bytes']
        values=[];assert len(expected[cid])==97
        for hour in range(97):
            r=read_cells(out/f'output{hour:08d}.xml',registry)
            assert r==expected[cid][hour]and r['time_h']==hour,(cid,hour)
            assert r['live_targets']==sum(r['live_targets_by_type'].values())
            ag0=sum(n for name,n in r['live_targets_by_type'].items()if byname[name]['antigen']=='Ag0')
            agpos=r['live_targets']-ag0
            row=dict(r,case=cid,ratio=rec['ratio'],arm=rec['arm'],seed=rec['seed'],target_fraction_of_initial=r['live_targets']/rec['targets'],effector_fraction_of_initial=r['live_effectors']/rec['effectors'],living_antigen_negative_targets=ag0,living_antigen_positive_targets=agpos,antigen_positive_fraction_of_living_targets=agpos/r['live_targets']if r['live_targets']else None)
            hourly.append(row);values.append(row)
        assert values[0]['live_targets']==rec['targets']and values[0]['live_effectors']==rec['effectors']
        area=integral([r['live_targets']for r in values]);trajectories[(rec['ratio'],rec['arm'],rec['seed'])]=values
        metrics.append({'case':cid,'ratio':rec['ratio'],'arm':rec['arm'],'seed':rec['seed'],'targets_96h':values[-1]['live_targets'],'target_AUC_cell_h':area,'time_average_targets_0_96h':area/96,'normalized_target_AUC_hours':area/rec['targets'],'target_fraction_96h':values[-1]['target_fraction_of_initial'],'antigen_positive_targets_96h':values[-1]['living_antigen_positive_targets'],'antigen_positive_fraction_96h':values[-1]['antigen_positive_fraction_of_living_targets']})
        balances.append(verify_balance(rec,out,cfg/'settings.xml'))
        print(json.dumps({'case':cid,'numbered_states_verified':97}),flush=True)
    primary=[]
    for ratio in [2,5]:
        for arm in ['2d','3d','local_entry']:
            r=[x for x in metrics if x['ratio']==ratio and x['arm']==arm]
            primary.append({'ratio':ratio,'arm':arm,'targets_96h':stats([x['targets_96h']for x in r]),'time_average_targets_0_96h':stats([x['time_average_targets_0_96h']for x in r]),'target_fraction_96h':math.fsum(x['target_fraction_96h']for x in r)/3,'normalized_target_AUC_hours':math.fsum(x['normalized_target_AUC_hours']for x in r)/3})
    close_nested(primary,read(HERE/'EXPECTED_PRIMARY_COMPARISON.json'))
    contrasts=[]
    for ratio in [2,5]:
        for left,right in [('3d','2d'),('local_entry','2d'),('local_entry','3d')]:
            for seed in range(3):
                a=trajectories[ratio,left,seed];b=trajectories[ratio,right,seed];dif=[x['live_targets']-y['live_targets']for x,y in zip(a,b)]
                contrasts.append({'ratio':ratio,'left':left,'right':right,'seed':seed,'hourly_live_target_difference':dif,'time_average_difference_cells':integral(dif)/96,'normalized_AUC_difference_hours':integral(dif)/1200,'endpoint_difference_cells':dif[-1],'positive_saved_hours':[i for i,x in enumerate(dif)if x>0],'zero_saved_hours':[i for i,x in enumerate(dif)if x==0],'negative_saved_hours':[i for i,x in enumerate(dif)if x<0]})
    write(args.output/'HOURLY_TRAJECTORIES.json',hourly);write(args.output/'CASE_METRICS.json',metrics);write(args.output/'PRIMARY_COMPARISON.json',primary);write(args.output/'HOURLY_CONTRASTS.json',contrasts);write(args.output/'MEDIUM_BALANCE_VERIFICATION.json',balances)
    write(args.output/'REPLAY_VERIFICATION.json',{'passed':True,'central_cases':18,'numbered_cell_states':1746,'field_states':1782,'balance_rows':3528,'raw_files_hash_verified':rawfiles,'raw_bytes_hash_verified':rawbytes,'portable_input_files_hash_verified':54,'scenario_ratio_groups':6,'seed_paired_contrast_trajectories':18,'previous_primary_comparison_reproduced':True,'engine_launched':False,'spatial_analysis_included':False,'scope':'Saved-data replay of central counts, assigned antigen composition, normalized hourly trajectories, integrals and IFNG/TNF medium inventories. Numerical convergence, biological calibration and final claim acceptance are evaluated separately.'})
    write(args.output/'ENVIRONMENT.json',{'python':sys.version,'implementation':platform.python_implementation(),'machine':platform.machine(),'byteorder':sys.byteorder,'dependencies':'Python standard library only'})
    write(args.output/'RESULT_MANIFEST_SHA256.json',{p.name:digest(p)for p in sorted(args.output.iterdir())if p.is_file()})
    print(json.dumps({'passed':True,'cases':18,'hourly_cell_states':1746,'field_states':1782}),flush=True)

if __name__=='__main__':main()
