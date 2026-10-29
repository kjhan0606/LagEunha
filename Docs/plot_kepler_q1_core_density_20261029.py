"""Common-scale native W0 density and initial-ID density change in the core."""
import argparse
import json
from pathlib import Path
import socket
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from plot_kepler_q1_boundary_comparison_20261029 import digest


def plot(initial,endpoints,output,summary):
    if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):
        raise RuntimeError('Render on a compute host')
    if output.exists() or summary.exists():raise FileExistsError('Refusing to overwrite a core plot')
    p0,m0=ordered_global(initial);ids=p0['u4if']&ID_MASK
    if m0['t']!=0:raise ValueError('Expected original IC')
    states=[(initial,p0,m0)]
    for path in endpoints:
        p,m=ordered_global(path,ids)
        if m['t']<=states[-1][2]['t'] or not np.array_equal(p['mass'],p0['mass']):
            raise ValueError('Expected increasing native checkpoints and fixed masses')
        states.append((path,p,m))
    axis=np.linspace(-.16,.16,650);x,y=np.meshgrid(axis,axis)
    xy=np.column_stack((x.ravel(),y.ravel()))
    fig,axes=plt.subplots(2,len(states),figsize=(5*len(states),9),squeeze=False,constrained_layout=True)
    records=[]
    for col,(path,p,meta) in enumerate(states):
        if np.any(p['w2']!=0) or not np.isfinite(p['den']).all() or np.any(p['den']<=0):
            raise ValueError('Requires finite positive W0 density')
        nearest=cKDTree(np.column_stack((p['x']-3.,p['y']-3.))).query(xy,workers=1)[1]
        fields=[p['den'][nearest].reshape(x.shape),100*(p['den'].astype(float)/p0['den'].astype(float)-1)[nearest].reshape(x.shape)]
        for row,field in enumerate(fields):
            ax=axes[row,col];norm=LogNorm(.004,.2) if row==0 else Normalize(-100,100)
            ax.imshow(field,origin='lower',extent=(-.16,.16,-.16,.16),norm=norm,cmap='viridis' if row==0 else 'coolwarm')
            ax.add_patch(Circle((0,0),.1,fill=False,color='white' if row==0 else 'black',ls='--',lw=.8))
            ax.set(title=f't = {meta["t"]:g}',xlabel='x - center',ylabel='y - center',aspect='equal')
        r=np.hypot(p['x']-3.,p['y']-3.);selected=p['den'][r<.15].astype(float)
        records.append({'time':meta['t'],'source':str(path),'sha256':digest(path),'current_r_lt_0p15_cells':len(selected),
                        'density_cell_percentiles_0_1_50_99_100':np.percentile(selected,[0,1,50,99,100]).tolist()})
    for row,label in enumerate(['Surface density (log scale)','Initial-ID density change (%)']):
        fig.colorbar(axes[row,-1].images[0],ax=axes[row,:].tolist(),label=label,extend='both',shrink=.9)
    fig.suptitle('Native bounded Q1 core: common scales, piecewise-constant W0 density')
    output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=170);plt.close(fig)
    result={'density_log_limits':[.004,.2],'density_change_percent_limits':[-100,100],
            'extent_centered':[-.16,.16,-.16,.16],'states':records,
            'scope':'Nearest-generator unsmoothed cell density. Dashed circle is current r=.1, while ID-relative change follows each generator from its own original density. Shared limits; colors saturate outside them. This figure alone cannot distinguish temporal, pressure, and boundary causes.'}
    with summary.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial',type=Path,required=True);parser.add_argument('--endpoints',type=Path,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--summary',type=Path,required=True)
    print(json.dumps(plot(**vars(parser.parse_args())),indent=2))
