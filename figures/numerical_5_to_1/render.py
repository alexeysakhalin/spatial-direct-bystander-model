"""Portable rendering of accepted TE5 numerical-followup data; no simulation."""
from pathlib import Path
import argparse,datetime,hashlib,json,math,os,resource
os.environ['MPLBACKEND']='Agg'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from PIL import Image

COLORS={'2d':'#D65114','local_entry':'#8F49C7'}
GRAY='#626973';LS=['-','--',':']
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def statistics(series):
 assert len(series)==3 and all(len(s)==97 for s in series)
 cols=list(zip(*series));return [[sum(v)/3 for v in cols],[min(v)for v in cols],[max(v)for v in cols]]
def panel(ax,group,col,letter,which=None):
 x=group['hours'];color=COLORS[group['arm']]
 titles=['Live targets','Paired target difference','Live effectors','Maximum field at live targets']
 labels=['Live targets (% of initial)','Refined - original (cells)','Live effectors (cells)','IFNG / TNF\n(model-relative level)']
 ax.set_title(letter+'  '+titles[col],loc='left',fontweight='bold',fontsize=10,pad=9)
 if col==1:
  ax.axhline(0,color='#9AA0A7',lw=.8,zorder=0)
  for seed in range(3):
   if which is not None and which!=f'seed{seed}':continue
   y=[a-b for a,b in zip(group['refined'][seed]['live_targets'],group['original'][seed]['live_targets'])]
   ax.plot(x,y,color=color,lw=1.35,ls=LS[seed],label=f'seed {seed}')
  ax.set_ylim(-90,90);ax.set_yticks([-80,-40,0,40,80])
 else:
  key={0:'live_targets',2:'live_effectors',3:'maximum_IFNG_TNF_at_living_targets'}[col]
  for role in ['original','refined']:
   if which is not None and which!=role:continue
   ys=[r[key]for r in group[role]]
   if col==0:ys=[[100*v/1200 for v in y]for y in ys]
   mean,lo,hi=statistics(ys);c=GRAY if role=='original'else color
   ax.fill_between(x,lo,hi,color=c,alpha=.13 if role=='original'else.18,lw=0)
   ax.plot(x,mean,color=c,ls='--'if role=='original'else'-',lw=1.8,label=role)
  if col==0:ax.set_ylim(0,103);ax.set_yticks([0,25,50,75,100])
  elif col==2:ax.set_ylim(200,520);ax.set_yticks([240,320,400,480])
  else:ax.set_ylim(0,group['field_plot_max']);ax.yaxis.set_major_locator(MaxNLocator(5))
 ax.set_xlim(0,96);ax.set_xticks([0,24,48,72,96]);ax.set_xlabel('Time (h)',fontsize=9);ax.set_ylabel(labels[col],fontsize=9,labelpad=5)
 ax.spines[['top','right']].set_visible(False);ax.spines[['bottom','left']].set_color('#9AA0A7');ax.grid(axis='y',color='#E9EBEF',lw=.55);ax.set_axisbelow(True);ax.tick_params(labelsize=8.5,length=3)
 ax.set_facecolor('white')
def render(input_dir,out):
 spec=load(input_dir/'PLOT_DATA_EN.json');provenance=load(input_dir/'PROVENANCE_EN.json')
 assert spec['seeds']==[0,1,2]and spec['initial_targets']==1200 and spec['initial_effectors']==240 and len(spec['groups'])==2
 for g in spec['groups']:
  assert g['hours']==list(range(97))
  for role in ['original','refined']:
   assert [r['seed']for r in g[role]]==[0,1,2]
   for s in g[role]:
    assert s['live_targets'][0]==1200 and s['live_effectors'][0]==240
    for k in ['live_targets','live_effectors','maximum_IFNG_TNF_at_living_targets']:assert len(s[k])==97 and all(math.isfinite(v)and v>=0 for v in s[k])
    assert all(g['field_plot_max']>=v for v in s['maximum_IFNG_TNF_at_living_targets'])
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42,'axes.unicode_minus':False})
 assert not out.exists();out.mkdir(parents=True);(out/'panels').mkdir();(out/'curves').mkdir();(out/'legends').mkdir()
 metadata={'CreationDate':datetime.datetime.fromisoformat(provenance['created_utc']),'ModDate':datetime.datetime.fromisoformat(provenance['created_utc']),'Creator':None,'Producer':None,'Title':'TE5 numerical followups: full 96 h trajectories'}
 fig,axs=plt.subplots(2,4,figsize=(14.4,8.6));fig.subplots_adjust(left=.065,right=.982,bottom=.17,top=.80,wspace=.33,hspace=.57)
 fig.suptitle('TE5 numerical followups: full 96 h trajectories',x=.065,y=.975,ha='left',fontsize=17,fontweight='bold')
 fig.text(.065,.927,'Initial targets:effectors = 5:1 (1,200:240)  |  Three matched seeds: 0, 1, 2',fontsize=11,color='#4C535D')
 handles=[Line2D([0],[0],color=GRAY,lw=1.8,ls='--',label='Original setting'),Line2D([0],[0],color='#343A40',lw=1.8,label='Refined setting (row color)')]+[Line2D([0],[0],color='#343A40',lw=1.35,ls=LS[i],label=f'Paired difference: seed {i}')for i in range(3)]
 fig.legend(handles=handles,ncol=5,frameon=False,bbox_to_anchor=(.06,.89),loc='lower left',fontsize=9.2,handlelength=2.5,columnspacing=1.4)
 letters='ABCDEFGH'
 for row,g in enumerate(spec['groups']):
  fig.text(.065,.854 if row==0 else.468,g['title'],fontsize=11.5,fontweight='bold',color=COLORS[g['arm']])
  for col in range(4):panel(axs[row,col],g,col,letters[row*4+col])
 fig.text(.065,.10,'Lines: three-seed means; shaded bands: min-max, not confidence intervals. B and F show same-seed differences.',fontsize=9.5)
 fig.text(.065,.074,'IFNG and TNF are identical model-relative fields here. The maximum is sampled at currently living targets.',fontsize=9.5)
 fig.text(.065,.048,'Completed numerical checks; whole-model convergence and biological calibration are not established.',fontsize=9.5,fontweight='bold',color='#49515A')
 fig.savefig(out/'TE5_NUMERICAL_FOLLOWUPS_N3.pdf',metadata=metadata);fig.savefig(out/'TE5_NUMERICAL_FOLLOWUPS_N3_600dpi.png',dpi=600);plt.close(fig)
 for row,g in enumerate(spec['groups']):
  for col in range(4):
   letter=letters[row*4+col];fig,ax=plt.subplots(figsize=(5.2,3.5));fig.subplots_adjust(left=.18,right=.97,bottom=.18,top=.82);panel(ax,g,col,letter);fig.suptitle(g['title'],x=.18,y=.98,ha='left',fontsize=10,color=COLORS[g['arm']]);fig.savefig(out/'panels'/f'{letter}_{g["arm"]}_600dpi.png',dpi=600);plt.close(fig)
   choices=[f'seed{i}'for i in range(3)]if col==1 else['original','refined']
   for choice in choices:
    fig,ax=plt.subplots(figsize=(5.2,3.5));fig.subplots_adjust(left=.18,right=.97,bottom=.18,top=.82);panel(ax,g,col,letter,choice);fig.suptitle(g['title']+' | '+choice,x=.18,y=.98,ha='left',fontsize=10,color=COLORS[g['arm']]);fig.savefig(out/'curves'/f'{letter}_{choice}_600dpi.png',dpi=600);plt.close(fig)
 for name,hs,width in [('setting',handles[:2],6.6),('seeds',handles[2:],8.4)]:
  fig=plt.figure(figsize=(width,.65));fig.legend(handles=hs,loc='center',ncol=len(hs),frameon=False,fontsize=10);fig.savefig(out/'legends'/f'{name}_600dpi.png',dpi=600);plt.close(fig)
 images=[]
 for f in sorted(out.rglob('*.png')):
  with Image.open(f)as im:
   assert all(abs(v-600)<.1 for v in im.info['dpi']);assert im.size[0]>=3000 or f.parent.name=='legends';images.append({'path':str(f.relative_to(out)),'size_px':list(im.size),'dpi':list(im.info['dpi']),'sha256':sha(f)})
 assert len(images)==29,len(images)
 save(out/'EXPORT_REVIEW_EN.json',{'render_complete':True,'PNG_count':29,'panel_PNG_count':8,'curve_PNG_count':18,'legend_PNG_count':2,'all_PNG_600dpi':True,'PDF_count':1,'images':images,'actual_visual_review_complete':False,'scientific_convergence_accepted':False})
 save(out/'MANIFEST_SHA256_EN.json',{str(f.relative_to(out)):sha(f)for f in sorted(out.rglob('*'))if f.is_file()})
 print(json.dumps({'out':str(out),'PNG_600dpi':29,'panels':8,'individual_curve_assets':18,'PDF':1,'convergence_accepted':False}))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:assert os.environ.get(k)=='1'
 resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(3*2**30,)*2);resource.setrlimit(resource.RLIMIT_CPU,(240,260));render(a.input_dir,a.out)
