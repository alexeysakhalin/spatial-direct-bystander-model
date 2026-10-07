from pathlib import Path
import json,xml.etree.ElementTree as ET
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from saved_frames import read_frame
CYTOKINES=("IFN-gamma","TNF","CXCL9-11","IL-2")

def initial_audit(out,rec,registry):
    time,d,types=read_frame(Path(out)/'initial.xml',registry)
    expected=np.loadtxt(Path(rec['folder'])/'cells.csv',delimiter=',',ndmin=2)
    assert expected.ndim==2 and expected.shape[1]==4 and np.isfinite(expected).all(),'invalid initial CSV'
    assert abs(time)<1e-10 and (d.dead==0).all(),'initial clock/death state'
    assert len(d)==len(expected),'initial cell count differs from exact CSV'
    actual=d[['position_0','position_1','position_2','cell_type']].to_numpy(dtype=float)
    assert (expected[:,3]==np.floor(expected[:,3])).all(),'noninteger CSV type'
    # Compare multisets; row order and ID allocation do not define cell identity.
    def ordered(a):return a[np.lexsort((a[:,3],a[:,2],a[:,1],a[:,0]))]
    assert np.allclose(ordered(actual),ordered(expected),rtol=0,atol=1e-7),'initial positions/types differ from CSV'
    reg={s['id']:s for s in json.loads(Path(registry).read_text())['states']}
    role=d.cell_type.map({k:s['role']for k,s in reg.items()})
    counts={k:int((role==k).sum()) for k in ('target','effector')}
    assert counts=={'target':rec['targets'],'effector':rec['effectors']},'declared initial role count mismatch'
    return {'passed':True,'initial_cells':len(d),'role_counts':counts,'exact_CSV_positions_and_types_matched':True,'time_h':time}

def audit_frame(xml,registry,no_cart=False):
    xml=Path(xml);t,d,types=read_frame(xml,registry)
    root=ET.parse(xml);variables=root.findall('.//microenvironment/domain/variables/variable')
    names={x.get('name'):4+int(x.get('ID'))for x in variables}
    assert len(names)==len(variables) and len(set(names.values()))==len(variables),'duplicate variable metadata'
    assert set(names)==set(CYTOKINES)|{'oxygen','debris'},'unexpected variable set'
    arrays=[v for k,v in loadmat(xml.with_name(xml.stem+'_microenvironment0.mat')).items()if not k.startswith('__')]
    assert len(arrays)==1,'ambiguous field matrix';a=arrays[0]
    assert a.ndim==2 and a.shape[0]==4+len(names) and np.isfinite(a).all(),'invalid field matrix'
    vol=a[3];assert (vol>0).all(),'invalid voxel volume'
    reg={s['id']:s for s in json.loads(Path(registry).read_text())['states']}
    target=d.cell_type.map({k:s['role']=='target'for k,s in reg.items()}).to_numpy(dtype=bool)
    live=(d.dead==0).to_numpy(dtype=bool)
    if no_cart:
        assert target.all(),'effector appeared in no-CAR control, including dead records'
        assert not any(reg[int(k)].get('primed',False)for k in d.cell_type),'priming appeared without CAR/IFNg source'
    points=d.loc[live&target,['position_0','position_1','position_2']].to_numpy()
    nearest=cKDTree(a[:3].T).query(points,k=1)[1] if len(points) else np.array([],dtype=int)
    row={'snapshot':xml.name,'time_h':t,'live_targets':int((live&target).sum()),'live_effectors':int((live&~target).sum()),'retained_dead':int((~live).sum()),'fields':{}}
    for name,index in names.items():
        c=a[index]
        assert c.min()>=-1e-12,(name,'negative concentration')
        if name=='oxygen':assert np.max(np.abs(c-38))<1e-7,'oxygen invariant'
        else:assert c.max()<=1+1e-8,(name,'above source-saturation bound')
        if no_cart and name in CYTOKINES:assert np.max(np.abs(c))<=1e-12,(name,'induced cytokine in no-CAR control')
        local=c[nearest]
        row['fields'][name]={'domain_volume_weighted_mean':float(np.dot(c,vol)/vol.sum()),'domain_min':float(c.min()),'domain_max':float(c.max()),'target_nearest_voxel_mean':float(local.mean())if len(local)else None,'target_nearest_voxel_quantiles_0_25_50_75_100':np.quantile(local,[0,.25,.5,.75,1]).tolist()if len(local)else None,'sampled_live_targets':len(local)}
    return row
