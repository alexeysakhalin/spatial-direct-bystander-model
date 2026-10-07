"""Figure7 print layout with three accepted realizations and preselected seed0 snapshots.
Verify-only mode creates no figure. Missing realizations cannot produce N3 plots.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,resource,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
import make_print_assets_v2 as base
from read_registered_frame import read_frame
F=base.F;ARMS=base.ARMS;C=base.C;HOURS=base.HOURS
COUNT_MAP={'total':'live_targets','positive':'live_Ag_positive','negative':'live_Ag_negative','effector':'live_effectors'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def verify_group(g,cases,ratio,arm):
 assert g['passed'] and g['n']==3 and g['seeds']==[0,1,2] and g['ratio']==ratio and g['arm']==arm
 assert len(cases)==3 and len(g['rows'])==97,'Exactly three complete realizations required'
 data={}
 for seed,c in enumerate(cases):
  assert c['passed'] and c['seed']==seed and c['ratio']==ratio and c['arm']==arm and c['id']==g['case_ids'][seed]
  a=c['acceptance_record'];assert a['exit_code']==0 and a['frames']==97 and a['medium_audit_passed']
  assert len(c['rows'])==97 and all(abs(r['time_h']-h)<1e-7 for h,r in enumerate(c['rows']))
  assert c['rows'][0]['live_targets']==1200 and c['rows'][0]['live_effectors']==1200//ratio
  assert all(r['live_Ag_positive']+r['live_Ag_negative']==r['live_targets']for r in c['rows'])
 for name,key in COUNT_MAP.items():
  y=np.array([[r[key]for r in c['rows']]for c in cases],dtype=float);assert np.isfinite(y).all()and(y>=0).all()and(y==np.floor(y)).all()
  mean=y.mean(axis=0);low=y.min(axis=0);high=y.max(axis=0)
  for h,r in enumerate(g['rows']):
   assert abs(r['time_h']-h)<1e-7;s=r['counts'][key];assert s['n_available']==3
   assert abs(s['mean']-mean[h])<1e-10 and s['min']==low[h] and s['max']==high[h],'Stored ensemble summary disagrees with individual runs'
  data[name]={'individual':y,'mean':mean,'min':low,'max':high}
 return data

def load_data(ratio,ensemble,input_index,registry):
 root=F.parent
 for name,digest in read(F/'snapshot_checksums.json').items():assert sha(root/name)==digest
 plots=read(F/f'counts_{ratio}to1.json');reg={s['id']:s for s in read(registry)['states']};data={};sources={}
 for arm in ARMS:
  v=plots[arm];curves={k:{n:np.asarray(z,dtype=float)for n,z in d.items()}for k,d in v['curves'].items()}
  for values in curves.values():
   y=values['individual'];assert y.shape==(3,97)
   assert np.array_equal(y.mean(0),values['mean'])and np.array_equal(y.min(0),values['min'])and np.array_equal(y.max(0),values['max'])
  frames={}
  for h in HOURS:
   xml=F/'snapshots'/str(ratio)/arm/f'output{h:08d}.xml'
   t,d,_=read_frame(xml,registry);assert abs(t-h)<1e-7
   target=d.cell_type.map({i:s['role']=='target'for i,s in reg.items()});live=d.dead==0
   masks={'positive':live&d.cell_type.map({i:s.get('antigen')=='AgHi'for i,s in reg.items()}),'negative':live&d.cell_type.map({i:s.get('antigen')=='Ag0'for i,s in reg.items()}),'effector':live&~target,'dead':~live}
   assert int((target&live).sum())==curves['total']['individual'][0,h]
   assert all(int(masks[k].sum())==curves[k]['individual'][0,h]for k in ['positive','negative','effector'])
   frames[h]=(d,masks)
  data[arm]={'curves':curves,'frames':frames,'seed0_id':v['case_ids'][0],'case_ids':v['case_ids']}
 return data,sources

def curves(ax,data,ymax):
 t=np.arange(97)
 for key in['total','positive','negative','effector']:
  d=data[key];style='--'if key=='effector'else'-'
  ax.fill_between(t,d['min'],d['max'],color=C[key],alpha=.10,lw=0)
  for run in d['individual']:ax.plot(t,run,color=C[key],ls=style,lw=.65,alpha=.40)
  ax.plot(t,d['mean'],color=C[key],ls=style,lw=1.9 if key=='total'else 1.6)
 ax.set(xlim=(0,96),ylim=(0,ymax),xticks=[0,24,48,72,96],yticks=np.arange(0,ymax+1,600),xlabel='Time (h)',ylabel='Living cells');ax.tick_params(labelsize=9);ax.grid(axis='y',color='#DCE3EC',lw=.65)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--ratio',type=int,choices=[2,5],required=True);ap.add_argument('--out',type=Path);ap.add_argument('--verify-only',action='store_true');ap.add_argument('--ensemble',type=Path,default=F/'ensemble_analysis_v1');ap.add_argument('--input-index',type=Path);ap.add_argument('--registry',type=Path,default=F/'state_registry.json');args=ap.parse_args()
 index=args.input_index or F/('TE2_EXISTING_CASES_EN.json'if args.ratio==2 else'MODEL_INPUT_INDEX_EN.json')
 # All acceptance, completeness and source checks precede output creation.
 data,sources=load_data(args.ratio,args.ensemble,index,args.registry)
 if args.verify_only:print(json.dumps({'passed':True,'ratio':args.ratio,'realizations':9,'groups':3,'numbered_states':873,'figure_created':False}));return
 assert args.out is not None and not args.out.exists();out=args.out;out.mkdir(parents=True);assets=out/'individual_png_600dpi';assets.mkdir()
 resource.setrlimit(resource.RLIMIT_CPU,(180,190));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
 ymax=max(2400,600*np.ceil(max(float(v['individual'].max())for d in data.values()for v in d['curves'].values())/600));exports=[];tag=f'TE{args.ratio}to1_N3'
 def export(fig,p,kind):
  fig.savefig(p,dpi=600,facecolor='white',metadata={'Software':None});plt.close(fig)
  with Image.open(p)as im:exports.append(dict(file=str(p.relative_to(out)),kind=kind,width_px=im.width,height_px=im.height,dpi=list(im.info.get('dpi',[])),sha256=sha(p)))
 for arm,letter in zip(ARMS,base.LABELS):
  for h in HOURS:
   fig=plt.figure(figsize=(3.2,3.2));ax=fig.add_axes([.055,.055,.89,.89],projection='3d'if arm=='3d'else None);base.snapshot(ax,arm,h,data[arm]['frames'][h],scale=1.45);export(fig,assets/f'{letter}_{arm}_{h:03d}h_TE{args.ratio}to1_SEED0.png','snapshot')
  fig,ax=plt.subplots(figsize=(3.7,3.2));fig.subplots_adjust(left=.22,right=.97,bottom=.19,top=.96);curves(ax,data[arm]['curves'],ymax);export(fig,assets/f'{letter}_{arm}_counts_{tag}.png','ensemble_curve')
 cellhandles=[Line2D([],[],ls='',marker='o',markersize=7,color=C[k],label=lab)for k,lab in [('positive','Ag+ target'),('negative','Ag- target'),('effector','CAR-T'),('dead','Retained dead cell')]]
 linehandles=[Line2D([],[],color=C[k],ls='--'if k=='effector'else'-',lw=2,label=lab)for k,lab in [('total','All targets'),('positive','Ag+ targets'),('negative','Ag- targets'),('effector','CAR-T')]]
 uncertainty=[Line2D([],[],color='#455367',lw=.7,alpha=.5,label='Individual run'),Line2D([],[],color='#455367',lw=2,label='Mean'),Patch(facecolor='#455367',alpha=.15,label='Min-max (n = 3)')]
 for name,handles in [('cell_legend',cellhandles),('curve_legend',linehandles),('realization_legend',uncertainty)]:
  fig=plt.figure(figsize=(10,.6));fig.legend(handles=handles,loc='center',ncol=len(handles),frameon=False);export(fig,assets/f'{name}_{tag}.png','legend')
 fig=plt.figure(figsize=(12,10.8),facecolor='white');fig.text(.045,.976,f'Targets:CAR-T = {args.ratio}:1 (1,200:{1200//args.ratio:,})',fontsize=15,weight='bold',va='top');fig.legend(handles=cellhandles,loc='upper right',bbox_to_anchor=(.978,.981),ncol=4,frameon=False,fontsize=10,columnspacing=1,handletextpad=.3)
 gs=fig.add_gridspec(3,5,left=.045,right=.980,top=.890,bottom=.140,width_ratios=[1,1,1,.24,1.18],wspace=.10,hspace=.35)
 for row,(arm,letter)in enumerate(zip(ARMS,base.LABELS)):
  anchor=gs[row,0].get_position(fig);fig.text(.045,anchor.y1+.018,letter+'   '+base.TITLES[arm],fontsize=14,weight='bold',va='bottom')
  for col,h in enumerate(HOURS):base.snapshot(fig.add_subplot(gs[row,col],projection='3d'if arm=='3d'else None),arm,h,data[arm]['frames'][h])
  curves(fig.add_subplot(gs[row,4]),data[arm]['curves'],ymax)
 fig.legend(handles=linehandles,loc='lower center',bbox_to_anchor=(.52,.071),ncol=4,frameon=False,fontsize=11)
 fig.legend(handles=uncertainty,loc='lower center',bbox_to_anchor=(.52,.032),ncol=3,frameon=False,fontsize=10)
 pdf=out/f'FIGURE7_{tag}.pdf';fig.savefig(pdf,dpi=600,metadata={'CreationDate':None,'ModDate':None,'Creator':None,'Producer':None});export(fig,out/f'FIGURE7_{tag}_600dpi.png','assembled_figure')

 plot_data={arm:{'case_ids':v['case_ids'],'snapshots_seed':0,'curves':{k:{n:a.tolist()for n,a in x.items()}for k,x in v['curves'].items()}}for arm,v in data.items()};save(out/'PLOT_DATA_EN.json',plot_data)
 report=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),ratio_targets_to_CAR_T=args.ratio,realizations_per_scenario=3,snapshot_seed=0,uncertainty='Pointwise descriptive min-max; not a confidence interval.',source_sha256=sources,builder_sha256=sha(__file__),snapshot_renderer_sha256=sha(Path(base.__file__)),exports=exports,pdf_sha256=sha(pdf),plot_data_sha256=sha(out/'PLOT_DATA_EN.json'),figure_render_completed=True)
 save(out/'PRINT_EXPORT_MANIFEST_EN.json',report);print(json.dumps({'output':str(out),'pdf':str(pdf),'png_files':len(exports),'visual_review_complete':False}))
if __name__=='__main__':main()
