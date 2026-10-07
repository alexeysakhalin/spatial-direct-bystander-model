from pathlib import Path
import hashlib,json,xml.etree.ElementTree as ET
import numpy as np
from scipy.io import loadmat
from saved_frames import read_frame

def load(p):return json.loads(Path(p).read_text())

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def quant(a):return np.quantile(a,[0,.25,.5,.75,1]).tolist() if len(a) else None

def reconstruct_gradient(c,xyz):
    """x-fast Cartesian storage; centered interior, one-sided first-order edges."""
    axes=[np.unique(a)for a in xyz]
    nx,ny,nz=map(len,axes);assert len(c)==nx*ny*nz
    z,y,x=np.meshgrid(axes[2],axes[1],axes[0],indexing='ij')
    assert np.array_equal(xyz,np.array([x.ravel(),y.ravel(),z.ravel()])), 'Cartesian ordering differs'
    cube=c.reshape(nz,ny,nx);g=[];spacing=[]
    for dim,axis in enumerate(axes):
        if len(axis)==1:
            g.append(np.zeros_like(c));spacing.append(None);continue
        delta=np.diff(axis);assert np.allclose(delta,delta[0],rtol=0,atol=1e-10)
        spacing.append(float(delta[0]))
        g.append(np.gradient(cube,float(delta[0]),axis=2-dim,edge_order=1).ravel())
    return np.array(g).T,axes,spacing

def voxel_indices(points,axes,spacing):
    ijk=[]
    for dim,(axis,dx)in enumerate(zip(axes,spacing)):
        if dx is None:ijk.append(np.zeros(len(points),dtype=int));continue
        lower=axis[0]-.5*dx;upper=axis[-1]+.5*dx
        # Actual reflective-domain outputs must be inside; do not emulate unsigned
        # underflow for positions that are not legitimate in these accepted runs.
        assert (points[:,dim]>=lower-1e-8).all()and(points[:,dim]<=upper+1e-8).all()
        ijk.append(np.clip(np.floor((points[:,dim]-lower)/dx).astype(int),0,len(axis)-1))
    i,j,k=ijk
    return (k*len(axes[1])+j)*len(axes[0])+i

def normalized(g):
    # BioFVM_vector.cpp normalize(std::vector<double>*): numerical regularizer,
    # not a biological receptor-occupancy or concentration-detection threshold.
    norm=np.sqrt(1e-32+g[:,0]**2+g[:,1]**2+g[:,2]**2)
    v=g/norm[:,None];v[norm<=1e-16]=0.
    return v

def read_transport(variable):
    assert variable is not None and variable.get('name')=='CXCL9-11'
    assert variable.get('units')=='dimensionless'
    dn=variable.find('./physical_parameter_set/diffusion_coefficient')
    kn=variable.find('./physical_parameter_set/decay_rate')
    assert dn is not None and kn is not None
    assert dn.get('units')=='micron^2/min' and kn.get('units')=='1/min'
    diffusion=float(dn.text);loss=float(kn.text)
    assert np.isfinite(diffusion) and np.isfinite(loss) and diffusion>=0 and loss>=0
    return diffusion,loss

def transport_scales(diffusion,loss):
    assert np.isfinite(diffusion) and np.isfinite(loss) and diffusion>=0 and loss>=0
    # With zero first-order removal there is no finite decay-only half-life or
    # diffusion-loss length. Null is retained; it is not a zero time/length.
    return (float(np.log(2)/loss),float(np.sqrt(diffusion/loss))) if loss>0 else (None,None)

def inspect(rec,accepted,registry):
    assert accepted['exit_code']==0 and accepted['frames']==97 and accepted['end_h']==96 and accepted['medium_audit_passed']
    for name,h in rec['input_sha256'].items():assert sha(Path(rec['folder'])/name)==h
    tree=ET.parse(Path(rec['folder'])/'settings.xml')
    assert tree.findtext('./microenvironment_setup/options/calculate_gradients')=='true'
    variable=tree.find('./microenvironment_setup/variable[@name="CXCL9-11"]')
    diffusion,loss=read_transport(variable)
    half_life,length=transport_scales(diffusion,loss)
    states={s['id']:s for s in load(registry)['states']};out=Path(rec['output']);rows=[];hashes={}
    for h in range(97):
        xml=out/f'output{h:08d}.xml';t,d,_=read_frame(xml,registry);assert t==h
        vars=ET.parse(xml).findall('.//microenvironment/domain/variables/variable')
        matches=[v for v in vars if v.get('name')=='CXCL9-11'];assert len(matches)==1
        field=xml.with_name(xml.stem+'_microenvironment0.mat')
        arrays=[v for k,v in loadmat(field).items()if not k.startswith('__')];assert len(arrays)==1
        a=arrays[0];assert a.shape[0]==4+len(vars) and np.isfinite(a).all()
        c=a[4+int(matches[0].get('ID'))];assert(c>=-1e-12).all()
        grad,axes,dx=reconstruct_gradient(c,a[:3]);norm=np.linalg.norm(grad,axis=1)
        target=d.cell_type.map({k:s['role']=='target'for k,s in states.items()}).to_numpy(bool)
        live=(d.dead==0).to_numpy(bool);em=live&~target;tm=live&target
        points=d.loc[em,['position_0','position_1','position_2']].to_numpy()
        indices=voxel_indices(points,axes,dx);gn=norm[indices]
        direction_norm=np.linalg.norm(normalized(grad[indices]),axis=1)
        counts={'live_targets':int(tm.sum()),'live_effectors':int(em.sum()),'retained_dead':int((~live).sum())}
        antigen=d.cell_type.map({k:s.get('antigen')for k,s in states.items()})
        primed=d.cell_type.map({k:s.get('primed',False)for k,s in states.items()}).to_numpy(bool)
        counts.update(live_Ag_positive=int((tm&antigen.isin(['AgHi','AgLo'])).sum()),live_Ag_negative=int((tm&(antigen=='Ag0')).sum()),live_primed_targets=int((tm&primed).sum()),live_unprimed_targets=int((tm&~primed).sum()))
        rows.append({'time_h':t,'counts':counts,'CXCL_domain_concentration_quantiles':quant(c),'CXCL_domain_gradient_magnitude_quantiles_per_um':quant(norm),'living_effector_local_concentration_quantiles':quant(c[indices]),'living_effector_local_gradient_magnitude_quantiles_per_um':quant(gn),'living_effector_direction_norm_quantiles':quant(direction_norm),'effectors_with_exact_zero_gradient':int((gn==0).sum()),'effectors_with_numerically_zero_direction':int((direction_norm==0).sum()),'effectors_with_direction_norm_at_least_0_99':int((direction_norm>=.99).sum()),'effectors_with_gradient_below_numeric_diagnostics':{str(cut):int((gn<cut).sum())for cut in [1e-16,1e-12,1e-8]}})
        for f in [xml,xml.with_name(xml.stem+'_cells.mat'),field]:hashes[str(f)]=sha(f)
    assert rows[-1]['counts']['live_targets']==accepted['last']['live_targets']
    return {'passed':True,'id':rec['id'],'arm':rec['arm'],'seed':rec['seed'],'rows':rows,'selected_times':{str(h):rows[h]for h in[0,1,6,24,48,72,96]},'grid_xyz':[len(x)for x in axes],'spacing_xyz_um':dx,'D_um2_min':diffusion,'loss_per_min':loss,'decay_only_half_life_min':half_life,'homogeneous_diffusion_loss_length_um':length,'input_sha256':rec['input_sha256'],'raw_sha256':hashes}

def selftests():
    checks=[];axes=[np.arange(4.)*2,np.arange(3.)*3,np.arange(3.)*4]
    z,y,x=np.meshgrid(axes[2],axes[1],axes[0],indexing='ij');xyz=np.array([x.ravel(),y.ravel(),z.ravel()])
    g,ax,dx=reconstruct_gradient(2*xyz[0]-3*xyz[1]+.5*xyz[2],xyz)
    assert np.allclose(g,[2,-3,.5]);checks.append('3D_linear_field_exact_interior_and_endcaps')
    assert np.array_equal(voxel_indices(xyz.T,ax,dx),np.arange(xyz.shape[1]));checks.append('all_voxel_centres_map_to_storage_indices')
    zero,*_=reconstruct_gradient(np.ones(xyz.shape[1]),xyz);assert not zero.any();checks.append('constant_field_zero_gradient')
    z,y,x=np.meshgrid([0.],axes[1],axes[0],indexing='ij');xyz2=np.array([x.ravel(),y.ravel(),z.ravel()]);g,_,_=reconstruct_gradient(xyz2[0]+xyz2[1],xyz2)
    assert np.allclose(g,[1,1,0]);checks.append('2D_gradient_z_is_zero')
    assert np.all(normalized(np.zeros((1,3)))==0);checks.append('zero_vector_stays_zero')
    norms=np.linalg.norm(normalized(np.array([[1e-4,0,0],[1e-12,0,0],[1e-18,0,0]])),axis=1)
    assert norms[0]>.999999 and norms[1]>.999999 and 0<norms[2]<.011;checks.append('normalization_amplitude_loss_and_numerical_regularization')
    try:reconstruct_gradient(np.ones(xyz.shape[1]),xyz[:,::-1])
    except AssertionError:checks.append('shuffled_storage_rejected')
    else:raise AssertionError('bad coordinate ordering accepted')
    return checks
