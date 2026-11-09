"""Publication figures from preserved cell averages and independent references."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
root=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'legend.fontsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'font.family':'serif'})
a=np.load(root/'profiles.npz')
labels=['Sod','Blast','Noh','Shu–Osher','Lax','Rarefaction','Collision','Static contact','Moving contact']
for field,col,j,ylab in [('density',3,0,r'$\rho$'),('pressure',5,2,r'$P$'),('velocity',4,1,r'$u$')]:
 fig,axes=plt.subplots(3,3,figsize=(7.2,7.2))
 fig.subplots_adjust(left=.09,right=.985,bottom=.065,top=.985,wspace=.34,hspace=.27)
 for k,(ax,label) in enumerate(zip(axes.flat,labels)):
  d=a[f'data{k}'];left=a[f'left{k}'];right=a[f'right{k}'];x=a[f'x{k}'];ref=a[f'ref{k}'];y=d[:,col];vref=ref[:,j]
  segments=np.stack([np.c_[left,y],np.c_[right,y]],axis=1)
  ax.add_collection(LineCollection(segments,colors='#0072B2',linewidths=.8,label='Cells'))
  ax.plot(x,vref,'k--',lw=.7,label='Reference')
  ax.plot(.5*(left+right),y,'.',color='#0072B2',ms=1.8)
  ax.set_xlim(left.min(),right.max());ax.set_xlabel(r'$x$',labelpad=0);ax.set_ylabel(ylab,labelpad=0)
  low=min(y.min(),vref.min());high=max(y.max(),vref.max());span=max(high-low,abs(high)*.05,.05)
  ax.set_ylim(low-.08*span,high+(.95 if k==0 else .40)*span)
  ax.text(.04,.97,label,transform=ax.transAxes,ha='left',va='top')
  ax.tick_params(pad=1);ax.locator_params(axis='both',nbins=3)
  if k==0:ax.legend(loc='upper left',bbox_to_anchor=(.01,.86),frameon=False,handlelength=1.3,borderpad=0,labelspacing=.1)
 fig.savefig(root.parent.parent/f'fig_material_1d_{field}.pdf')
 fig.savefig(root.parent.parent/f'fig_material_1d_{field}.png',dpi=160)
 plt.close(fig)
