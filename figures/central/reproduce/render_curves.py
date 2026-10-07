from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
C={'positive':'#006AFF','negative':'#FF7700','effector':'#00A465','total':'#202B3B'}
N={'positive':'Ag+ targets','negative':'Ag- targets','effector':'CAR-T','total':'All targets'}
ARM={'2d':'2D monolayer','3d':'3D spheroid','local_entry':'Local entry'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):assert not p.exists();p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--ratio',type=int,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists();a.output.mkdir(parents=True)
 data=json.loads(a.source.read_text());plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False});t=np.arange(97);exports=[]
 for arm in ['2d','3d','local_entry']:
  for key in ['total','positive','negative','effector']:
   d=data[arm]['curves'][key];runs=np.asarray(d['individual'],dtype=float);assert runs.shape==(3,97);assert np.array_equal(runs.mean(0),d['mean'])and np.array_equal(runs.min(0),d['min'])and np.array_equal(runs.max(0),d['max'])
   for kind in ['seed0','seed1','seed2','mean','ensemble']:
    fig,ax=plt.subplots(figsize=(3.7,3.2));fig.subplots_adjust(left=.22,right=.97,bottom=.19,top=.86);style='--'if key=='effector'else'-';curves=[]
    if kind=='ensemble':
     ax.fill_between(t,d['min'],d['max'],color=C[key],alpha=.10,lw=0)
     for run in runs:curves+=ax.plot(t,run,color=C[key],ls=style,lw=.65,alpha=.40)
     curves+=ax.plot(t,d['mean'],color=C[key],ls=style,lw=1.9 if key=='total'else 1.6);expected=list(runs)+[np.asarray(d['mean'])]
    elif kind=='mean':curves+=ax.plot(t,d['mean'],color=C[key],ls=style,lw=1.9 if key=='total'else 1.6);expected=[np.asarray(d['mean'])]
    else:
     i=int(kind[-1]);curves+=ax.plot(t,runs[i],color=C[key],ls=style,lw=1.6);expected=[runs[i]]
    for line,values in zip(curves,expected):assert np.array_equal(line.get_xdata(),t)and np.array_equal(line.get_ydata(),values)
    subtitle={'mean':'Mean of 3 seeds','ensemble':'Seeds, mean and min-max'}.get(kind,'Seed '+kind[-1]);fig.text(.22,.965,f'{a.ratio}:1 | {ARM[arm]}',fontsize=9,va='top');fig.text(.22,.917,N[key]+' | '+subtitle,fontsize=9,va='top')
    ax.set(xlim=(0,96),ylim=(0,2400),xticks=[0,24,48,72,96],yticks=np.arange(0,2401,600),xlabel='Time (h)',ylabel='Living cells');ax.tick_params(labelsize=9);ax.grid(axis='y',color='#DCE3EC',lw=.65)
    p=a.output/f'{arm}_{key}_{kind}.png';fig.savefig(p,dpi=600,facecolor='white',metadata={'Software':None});plt.close(fig)
    with Image.open(p)as im:
     assert im.size==(2220,1920)and all(abs(v-600)<.02 for v in im.info['dpi'])and not im.info.get('Software');assert im.convert('RGB').getpixel((0,0))==(255,255,255)
    exports.append({'file':p.name,'ratio':a.ratio,'arm':arm,'variable':key,'curve':kind,'seeds':[0,1,2]if kind in ['mean','ensemble']else[int(kind[-1])],'dimensions_px':[2220,1920],'dpi':600,'sha256':sha(p),'plotted_values_checked':len(expected)*97})
 save(a.output/'assets.json',{'initial_target_to_effector_ratio':f'{a.ratio}:1','case_ids':{arm:data[arm]['case_ids']for arm in data},'axis_limits':{'hours':[0,96],'living_cells':[0,2400]},'data':'plot_data.json','assets':exports})
 (a.output/'plot_data.json').write_bytes(a.source.read_bytes());save(a.output/'checksums_sha256.json',{p.name:sha(p)for p in sorted(a.output.iterdir())if p.is_file()});print(json.dumps({'ratio':a.ratio,'PNG':len(exports)}))
if __name__=='__main__':main()
