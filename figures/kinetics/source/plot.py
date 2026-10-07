from pathlib import Path
import argparse,hashlib,json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'pdf.fonttype':42,'savefig.facecolor':'white','figure.facecolor':'white','axes.facecolor':'white'})
TEAL='#007D80';GRAY='#5D6269';LIGHT='#B8BDC3'
def panel(ax,r,letter=None,layer='all'):
 t=r['time_h'];b=r['baseline_normalized'];lo=[min(x)for x in zip(*b)];hi=[max(x)for x in zip(*b)]
 if layer in ['all','range']:ax.fill_between(t,lo,hi,color=LIGHT,alpha=.35,lw=0,zorder=1)
 for seed in [1,2,0]:
  if layer=='all'or layer=='baseline_seed'+str(seed):ax.plot(t,b[seed],color=GRAY if seed==0 else LIGHT,lw=1.2 if seed==0 else .65,ls='--'if seed==0 else '-',zorder=2)
 if layer in ['all','challenge']:ax.plot(t,r['challenge_normalized'],color=TEAL,lw=1.6,zorder=3)
 ax.set(xlim=(0,96),ylim=(0,1.08),xticks=[0,24,48,72,96],yticks=[0,.25,.5,.75,1],xlabel='Time (h)',ylabel='Living targets / initial targets');ax.tick_params(length=3,width=.6);ax.grid(axis='y',color='#E5E7E9',lw=.5,zorder=0)
 if layer=='all':ax.set_title((letter+'  'if letter else'')+r['representation'],loc='left',fontweight='bold');ax.text(.98,.94,r['challenge_label'],ha='right',va='top',transform=ax.transAxes,fontsize=8,color='#3B3E43')
 else:ax.set_axis_off();ax.grid(False)
def legend(ax):
 ax.set_axis_off();handles=[Line2D([],[],color=TEAL,lw=1.6,label='Challenge, seed 0'),Line2D([],[],color=GRAY,lw=1.2,ls='--',label='Baseline, seed 0'),Line2D([],[],color=LIGHT,lw=.7,label='Baseline, seeds 1 and 2'),Patch(color=LIGHT,alpha=.35,label='Baseline N3 min-max')];ax.legend(handles=handles,loc='center',ncol=2,frameon=False,fontsize=9,columnspacing=2,handlelength=2.8)
def png(fig,p,title):fig.savefig(p,dpi=600,facecolor='white',metadata={'Software':None,'Title':title})
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output;assert not out.exists();out.mkdir(parents=True);d=json.loads((ROOT/'data.json').read_text());assert len(d['cases'])==9 and d['initial_targets']==1200 and d['initial_effectors']==600
 for i,r in enumerate(d['cases']):
  assert len(r['time_h'])==len(r['challenge_normalized'])==97 and len(r['baseline_normalized'])==3;label=chr(65+i);stem=r['figure_id'];fig,ax=plt.subplots(figsize=(3.55,2.8));fig.subplots_adjust(left=.18,right=.97,bottom=.19,top=.86);panel(ax,r,label);png(fig,out/(stem+'_panel.png'),'Kinetic comparison '+label);plt.close(fig)
  for layer in ['challenge','baseline_seed0','baseline_seed1','baseline_seed2','range']:
   fig,ax=plt.subplots(figsize=(3.55,2.8));fig.subplots_adjust(left=.18,right=.97,bottom=.19,top=.86);panel(ax,r,layer=layer);png(fig,out/(stem+'_'+layer+'.png'),'Kinetic comparison '+label+' '+layer);plt.close(fig)
 fig,ax=plt.subplots(figsize=(7.1,.72));fig.subplots_adjust(left=0,right=1,bottom=0,top=1);legend(ax);png(fig,out/'legend.png','Kinetic comparison legend');plt.close(fig)
 fig,axes=plt.subplots(3,3,figsize=(10.65,8.55));fig.subplots_adjust(left=.068,right=.99,bottom=.115,top=.96,wspace=.3,hspace=.46)
 for i,(ax,r)in enumerate(zip(axes.flat,d['cases'])):panel(ax,r,chr(65+i))
 ax=fig.add_axes([.12,.005,.76,.077]);legend(ax);png(fig,out/'kinetic_comparisons.png','Growth and background-loss comparisons');fig.savefig(out/'kinetic_comparisons.pdf',metadata={'Title':'Growth and background-loss comparisons','Creator':None,'Producer':None,'CreationDate':None,'ModDate':None});plt.close(fig)
 manifest={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(out.rglob('*'))if p.is_file()};(out/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n');print(json.dumps({'panels':9,'curve_layers':36,'range_layers':9,'PNG_files':56,'PDF_files':1,'dpi':600}))
if __name__=='__main__':main()
