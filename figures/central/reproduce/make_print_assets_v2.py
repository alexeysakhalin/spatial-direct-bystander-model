from pathlib import Path
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
F=Path(__file__).resolve().parent.parent/'source_data'
C={'positive':'#006AFF','negative':'#FF7700','effector':'#00A465','dead':'#B5BCC6','total':'#202B3B'}
ARMS=['2d','3d','local_entry'];LABELS=['A','B','C'];HOURS=[0,24,96]
TITLES={'2d':'2D monolayer','3d':'3D spheroid','local_entry':'Local-entry tissue scenario'}
def snapshot(ax,arm,h,frame,scale=1):
 d,m=frame
 if arm=='3d':
  visible=~((d.position_0>0)&(d.position_1<0));cuts=[[(0,-235,-235),(0,0,-235),(0,0,235),(0,-235,235)],[(0,0,-235),(235,0,-235),(235,0,235),(0,0,235)]]
  ax.add_collection3d(Poly3DCollection(cuts,facecolors='#E1E8F0',edgecolors='#A9B6C6',linewidths=.55,alpha=.14))
  for k,s in [('dead',7),('positive',12),('negative',12),('effector',9)]:
   mask=m[k]&visible;ax.scatter(d.loc[mask,'position_0'],d.loc[mask,'position_1'],d.loc[mask,'position_2'],c=C[k],s=s*scale,alpha=.55 if k=='dead'else 1,depthshade=False,edgecolors='none',rasterized=False)
  ax.view_init(elev=23,azim=-60);ax.set(xlim=(-250,250),ylim=(-250,250),zlim=(-250,250));ax.set_box_aspect((1,1,1));ax.set_axis_off()
  ax.plot([-210,-110],[-235,-235],[-235,-235],color=C['total'],lw=2.4);ax.text(-160,-245,-245,'100 µm',fontsize=9,ha='center',va='top')
  ax.text2D(.04,.96,f'{h} h',transform=ax.transAxes,fontsize=12,weight='bold',va='top')
 else:
  bound=1180 if arm=='2d'else 600;visible=np.ones(len(d),dtype=bool)if arm=='2d'else (np.abs(d.position_2)<=60).to_numpy()
  for k,s in [('dead',3),('positive',5.6 if arm=='2d'else 13),('negative',5.6 if arm=='2d'else 13),('effector',4.9 if arm=='2d'else 10)]:
   mask=m[k]&visible;ax.scatter(d.loc[mask,'position_0'],d.loc[mask,'position_1'],c=C[k],s=s*scale,edgecolors='none',alpha=.55 if k=='dead'else 1,rasterized=False)
  ax.set(xlim=(-bound,bound),ylim=(-bound,bound),aspect='equal',xticks=[],yticks=[]);ax.set_facecolor('white')
  for sp in ax.spines.values():sp.set_visible(False)
  ax.text(.035,.965,f'{h} h',transform=ax.transAxes,fontsize=12,weight='bold',va='top',bbox={'facecolor':'white','edgecolor':'none','pad':1.5})
  xx=.35*bound;yy=-.80*bound;ax.plot([xx,xx+200],[yy,yy],lw=2.2,color=C['total']);ax.text(xx+100,yy+.035*bound,'200 µm',fontsize=9,ha='center',va='bottom')
  if arm=='local_entry':ax.text(.965,.965,'|z| ≤ 60 µm',transform=ax.transAxes,fontsize=8.5,ha='right',va='top',color='#475569')
