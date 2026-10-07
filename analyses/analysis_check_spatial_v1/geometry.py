from pathlib import Path
import hashlib,json
import numpy as np
from scipy.spatial import ConvexHull,cKDTree
from saved_frames import read_frame
REGISTRY=None
POSITION=['position_0','position_1','position_2']
TOL_UM=1e-8
FIRST_KEYS=['living_effectors_inside_initial_target_hull','living_effectors_touching_living_target_spheres','living_effectors_with_saved_live_target_attack_link','living_effectors_with_positive_recent_activity_memory']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def load(p):return json.loads(Path(p).read_text())

def quant(a):return np.quantile(a,[0,.25,.5,.75,1]).tolist() if len(a) else None

def inside(points,equations):
 # Half-space equations are normalized by scipy; tolerance is in micrometres.
 if not len(points):return np.zeros(0,dtype=bool)
 return ((points@equations[:,:-1].T+equations[:,-1])<=TOL_UM).all(axis=1)

def sphere_contacts(ep,er,tp,tr):
 assert len(ep)==len(er) and len(tp)==len(tr)
 if not len(ep) or not len(tp):return np.zeros(len(ep),dtype=bool),None
 tree=cKDTree(tp);distance=tree.query(ep,k=1,workers=1)[0]
 candidates=tree.query_ball_point(ep,er+float(tr.max())+TOL_UM,workers=1)
 touching=np.zeros(len(ep),dtype=bool)
 for i,js in enumerate(candidates):
  if js:
   j=np.asarray(js,dtype=int);touching[i]=bool((np.linalg.norm(tp[j]-ep[i],axis=1)<=er[i]+tr[j]+TOL_UM).any())
 return touching,quant(distance)

def links(frame,effector_mask,role):
 ids=frame.ID.to_numpy();assert np.isfinite(ids).all() and (ids==np.floor(ids)).all()
 target_ids=frame.loc[role=='target','ID'].to_numpy(dtype=int);live_ids=frame.loc[(role=='target')&(frame.dead==0),'ID'].to_numpy(dtype=int)
 a=frame.loc[effector_mask,'attack_target'].to_numpy();assert np.isfinite(a).all() and (a==np.floor(a)).all() and (a>=-1).all()
 present=a[a>=0].astype(int);assert np.isin(present,target_ids).all(),'Saved effector attack link points outside recorded target population'
 return int(np.isin(present,live_ids).sum()),int(len(present)-np.isin(present,live_ids).sum())

def inspect(rec,accepted):
 assert accepted['exit_code']==0 and accepted['frames']==97 and accepted['medium_audit_passed']
 root=Path(rec['output']);assert all(sha(Path(rec['folder'])/n)==h for n,h in rec['input_sha256'].items())
 states={r['id']:r for r in load(REGISTRY)['states']};_,init,_=read_frame(root/'initial.xml',REGISTRY)
 roles={i:s['role']for i,s in states.items()};ir=init.cell_type.map(roles);ip=init.loc[ir=='target',POSITION].to_numpy();dimension=2 if rec['arm']=='2d' else 3
 assert len(ip)==1200 and len(init.loc[ir=='effector']) in[0,240,600]
 hull=ConvexHull(ip[:,:dimension]);equations=hull.equations.copy();rows=[];hashes={str(root/'initial.xml'):sha(root/'initial.xml'),str(root/'initial_cells.mat'):sha(root/'initial_cells.mat')}
 for h in range(97):
  xml=root/f'output{h:08d}.xml';t,d,_=read_frame(xml,REGISTRY);assert abs(t-h)<1e-7
  for column in ['ikm_recent_activity','attack_target']:assert column in d and np.isfinite(d[column]).all()
  role=d.cell_type.map(roles);live=d.dead==0;em=(role=='effector')&live;tm=(role=='target')&live
  assert (d.total_volume>0).all() and (d.ikm_recent_activity>=-1e-12).all()
  ep=d.loc[em,POSITION].to_numpy();tp=d.loc[tm,POSITION].to_numpy();er=(3*d.loc[em,'total_volume'].to_numpy()/(4*np.pi))**(1/3);tr=(3*d.loc[tm,'total_volume'].to_numpy()/(4*np.pi))**(1/3)
  inn=inside(ep[:,:dimension],equations);touch,distance=sphere_contacts(ep,er,tp,tr);live_links,dead_links=links(d,em,role)
  row=dict(time_h=t,living_targets=int(tm.sum()),living_effectors=int(em.sum()),living_effectors_inside_initial_target_hull=int(inn.sum()),living_effectors_touching_living_target_spheres=int(touch.sum()),living_effectors_with_saved_live_target_attack_link=live_links,living_effectors_with_saved_dead_target_attack_link=dead_links,living_effectors_with_positive_recent_activity_memory=int((d.loc[em,'ikm_recent_activity']>0).sum()),nearest_living_target_centre_distance_quantiles_um=distance)
  row['fraction_living_effectors_inside_initial_target_hull']=float(inn.mean())if len(inn)else None
  rows.append(row)
  for x in [xml,xml.with_name(xml.stem+'_cells.mat')]:hashes[str(x)]=sha(x)
 first={k:next((r['time_h']for r in rows if r[k]>0),None)for k in FIRST_KEYS}
 result=dict(passed=True,id=rec['id'],arm=rec['arm'],seed=rec['seed'],rows=rows,first_saved_positive_hours=first,selected_times={str(h):rows[h]for h in[0,1,6,24,48,72,96]},initial_target_envelope={'dimension':dimension,'definition':'Convex hull of initial target centres;2D uses xy,3D uses xyz. It is fixed over time and ignores target radii.','halfspace_equations':equations.tolist(),'numerical_boundary_tolerance_um':TOL_UM},input_sha256=rec['input_sha256'],raw_cell_input_sha256=hashes,scope='Current living effector population includes descendants. Interior membership is not vessel crossing. Sphere overlap is geometric proximity, not proof of antigen recognition. Saved attack links are instantaneous engine states, not counted kills. Recent activity can persist after contact. First saved positive time is not the true first event; intervals between hourly outputs are unobserved. No migration speed, lineage tree or death-mechanism fraction inferred.')
 return result

def selftest():
 checks=[]
 p=np.array([[0.,0.,0.],[5,0,0.]])
 touching,q=sphere_contacts(p,np.ones(2),np.array([[2.,0.,0.]]),np.ones(1));assert touching.tolist()==[True,False];checks.append('touching_and_distant_spheres')
 # A larger, more distant target may overlap although the nearest centre does not.
 touching,_=sphere_contacts(np.zeros((1,3)),np.array([1.]),np.array([[3.,0,0],[4.,0,0]]),np.array([1.,4.]));assert touching[0];checks.append('all_candidate_radii_not_only_nearest_centre')
 z,q=sphere_contacts(p,np.ones(2),np.empty((0,3)),np.array([]));assert not z.any()and q is None;checks.append('no_targets_distance_undefined')
 eq=ConvexHull(np.array([[0.,0.],[1,0],[1,1],[0,1]])).equations;assert inside(np.array([[.5,.5],[1,0],[1.1,0]]),eq).tolist()==[True,True,False];checks.append('fixed_hull_inside_edge_outside')
 import pandas as pd
 f=pd.DataFrame({'ID':[1,2,3],'dead':[0,0,1],'attack_target':[2,-1,-1]});role=pd.Series(['effector','target','target']);em=role=='effector';assert links(f,em,role)==(1,0)
 f.loc[0,'attack_target']=3;assert links(f,em,role)==(0,1);checks.append('living_and_dead_attack_target_distinguished')
 f.loc[0,'attack_target']=99
 try:links(f,em,role)
 except AssertionError:checks.append('missing_attack_target_rejected')
 else:raise AssertionError('Invalid saved link accepted')
 return checks
