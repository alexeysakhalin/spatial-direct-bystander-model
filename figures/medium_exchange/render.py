"""Rebuild the matched medium-exchange figures from the included plotted data."""
from pathlib import Path
import argparse, json, os, tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'spatial_model_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
args=parser.parse_args()
OUT=args.output.resolve()
if OUT.exists():
    parser.error('Use a fresh output directory')
data=json.loads((Path(__file__).resolve().parent/'PLOT_DATA_EN.json').read_text())
OUT.mkdir(parents=True)
(OUT/'individual_png_600dpi').mkdir()
arms=['2d','3d','local_entry']
titles=['2D culture','3D spheroid','Tissue-entry scenario']
colors={0:'#5944A3',10:'#008B87'}
styles={0:'-',10:'--'}
cases={}
for arm in arms:
    for flow in [0,10]:
        rows=data['series'][f'turnover{flow}_{arm}_seed0']
        assert len(rows)==97 and [r['time_h'] for r in rows]==list(range(97))
        assert all(r['IFNg_domain_mean']==r['TNF_domain_mean'] for r in rows)
        cases[arm,flow]={'rows':[dict(r,fields={'IFN-gamma':{'domain_volume_weighted_mean':r['IFNg_domain_mean']}}) for r in rows]}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
target_max=max(r['live_targets'] for c in cases.values() for r in c['rows'])
ag0_max=max(r['live_Ag_negative'] for c in cases.values() for r in c['rows'])
field_max=max(r['fields']['IFN-gamma']['domain_volume_weighted_mean'] for c in cases.values() for r in c['rows'])
limits=[(0,max(1250,((target_max//250)+1)*250)),(0,((ag0_max//100)+1)*100),(0,((int(field_max*10)//1)+1)/10)]
labels=['Living targets','Living Ag− targets','Mean IFNγ / TNF\n(relative concentration)']
def panel(ax,arm,row):
    for flow in [0,10]:
        d=cases[arm,flow];x=[r['time_h'] for r in d['rows']]
        y=[r['live_targets'] if row==0 else r['live_Ag_negative'] if row==1 else r['fields']['IFN-gamma']['domain_volume_weighted_mean'] for r in d['rows']]
        ax.plot(x,y,color=colors[flow],linestyle=styles[flow],lw=2.1)
    ax.set(xlim=(0,96),ylim=limits[row],xticks=[0,24,48,72,96])
    ax.grid(axis='y',alpha=.15,color='#64748B');ax.set_axisbelow(True);ax.tick_params(length=3)
    ax.set_xlabel('Time (h)');ax.set_ylabel(labels[row])
fig,axs=plt.subplots(3,3,figsize=(10.4,8.6))
fig.subplots_adjust(left=.095,right=.985,bottom=.075,top=.845,hspace=.52,wspace=.40)
for col,(arm,title) in enumerate(zip(arms,titles)):
    for row in range(3):
        ax=axs[row,col];panel(ax,arm,row)
        if col:ax.set_ylabel('')
        if row==0:ax.set_title(title,pad=12)
        ax.text(-.20,1.05,chr(ord('A')+3*row+col),transform=ax.transAxes,fontweight='bold',fontsize=13)
        pf,pa=plt.subplots(figsize=(3.5,2.8));panel(pa,arm,row);pa.set_title(title)
        pf.subplots_adjust(left=.24,right=.96,bottom=.23,top=.83)
        pf.savefig(OUT/'individual_png_600dpi'/f'{arm}_{["targets","Ag_negative","cytokines"][row]}_600dpi.png',dpi=600,facecolor='white',metadata={'Software':None});plt.close(pf)
fig.text(.5,.975,'External medium exchange | Targets:CAR-T = 2:1',ha='center',va='top',fontsize=14)
handles=[Line2D([0],[0],color=colors[f],linestyle=styles[f],lw=2.2,label=label) for f,label in [(0,'Closed medium'),(10,'Exchange: 10 total volumes / 96 h')]]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.94),ncol=2,frameon=False,fontsize=10)
fig.savefig(OUT/'PAIRED_MEDIUM_EXCHANGE_TE2_SEED0.pdf',metadata={'Title':'Paired medium-exchange model contrasts','Creator':None,'Producer':None,'CreationDate':None})
fig.savefig(OUT/'PAIRED_MEDIUM_EXCHANGE_TE2_SEED0_600dpi.png',dpi=600,facecolor='white',metadata={'Software':None});plt.close(fig)
fig,ax=plt.subplots(figsize=(8,.6));ax.axis('off')
ax.legend(handles=handles,loc='center',ncol=2,frameon=False,fontsize=11)
fig.subplots_adjust(0,0,1,1)
fig.savefig(OUT/'individual_png_600dpi/legend_600dpi.png',dpi=600,facecolor='white',metadata={'Software':None})
plt.close(fig)
print(OUT)
