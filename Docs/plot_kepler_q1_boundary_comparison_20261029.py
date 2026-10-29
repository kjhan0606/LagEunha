"""Same-time density fields and ID-matched radial errors for bounded/periodic Q1.

Piecewise-constant W0 Voronoi pixels, common fixed color limits, no smoothing.
Use only complete native checkpoints. Geometry/rendering runs on a compute host.
"""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE, mean


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sample(p, values, pixels=720, bounded=False, disk_mask=None):
    axis = np.linspace(-2.25,2.25,pixels)
    x,y = np.meshgrid(axis,axis)
    xy = np.column_stack((p['x']-3.,p['y']-3.))
    index = cKDTree(xy).query(np.column_stack((x.ravel(),y.ravel())),workers=1)[1]
    image = np.asarray(values)[index].reshape(x.shape).copy()
    if bounded:
        angle = np.arctan2(y,x)
        k = np.rint(angle/(2*np.pi/1024)-.5)
        normal = (k+.5)*2*np.pi/1024
        image[x*np.cos(normal)+y*np.sin(normal)>2.75*np.cos(np.pi/1024)] = np.nan
    if disk_mask is not None:
        image[~disk_mask[index].reshape(x.shape)] = np.nan
    return image


def radial(initial, final):
    r = np.hypot(initial['x']-3.,initial['y']-3.)
    delta = final['den'].astype(float)/initial['den'].astype(float)-1
    entropy = np.log(final['ie'].astype(float)/initial['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(final['volume'].astype(float)/initial['volume'].astype(float))
    mass = initial['mass'].astype(float);rows=[]
    for lo,hi in zip(np.linspace(.5,2.,31)[:-1],np.linspace(.5,2.,31)[1:]):
        mask = (r>=lo)&(r<hi)
        if np.any(mask):
            rows.append({'lo':float(lo),'hi':float(hi),'density_abs':mean(abs(delta[mask]),mass[mask]),
                         'entropy_signed':mean(entropy[mask],mass[mask]),'entropy_abs':mean(abs(entropy[mask]),mass[mask])})
    return delta,rows


def plot(initial, endpoint, periodic_initial, periodic_endpoint, output, summary):
    if output.exists() or summary.exists():
        raise FileExistsError('Refusing to overwrite a density comparison')
    a,m0 = ordered_global(initial);ids = a['u4if']&ID_MASK
    b,mt = ordered_global(endpoint,ids)
    c,mc = ordered_global(periodic_initial,np.arange(384**2,dtype=np.uint64))
    d,md = ordered_global(periodic_endpoint,np.arange(384**2,dtype=np.uint64))
    if m0['t']!=mc['t'] or m0['t']!=0 or mt['t']!=md['t'] or mt['t']<=0:
        raise ValueError('Requires exact matching initial/final physical times')
    if not np.array_equal(a['mass'],b['mass']) or not np.array_equal(c['mass'],d['mass']):
        raise ValueError('Fixed mass gate failed')
    if any(np.any(q['w2']!=0) for q in [a,b,c,d]):
        raise ValueError('Nearest-generator rendering is W0 only')
    disk = (np.hypot(a['x']-3.,a['y']-3.)>=.5)&(np.hypot(a['x']-3.,a['y']-3.)<2)
    original_disk = (np.hypot(c['x']-3.,c['y']-3.)>=.5)&(np.hypot(c['x']-3.,c['y']-3.)<2)
    if not np.array_equal(ids[disk],(c['u4if']&ID_MASK)[original_disk]):
        raise ValueError('Scientific disk ID sets differ')
    if any(not np.array_equal(a[k][disk],c[k][ids[disk].astype(int)]) for k in ['x','y','vx','vy','mass','ie']):
        raise ValueError('Scientific disk initial data differ')
    delta,rb = radial(a,b);old_delta,rp = radial(c,d)
    images = [sample(a,a['den'],bounded=True),sample(b,b['den'],bounded=True),sample(d,d['den']),
              sample(b,100*delta,bounded=True,disk_mask=disk),sample(d,100*old_delta,disk_mask=original_disk)]
    fig,axes = plt.subplots(2,3,figsize=(15,9.5),constrained_layout=True)
    titles = ['Common disk initial state',f'Bounded Q1, t={mt["t"]:g}',f'Periodic Q1, t={md["t"]:g}',
              'Bounded: ID-matched density change (%)','Periodic: ID-matched density change (%)']
    for i,(ax,field,title) in enumerate(zip(axes.flat,images,titles)):
        norm = Normalize(0,1.05) if i<3 else Normalize(-1,1)
        im = ax.imshow(field,origin='lower',extent=(-2.25,2.25,-2.25,2.25),norm=norm,cmap='viridis' if i<3 else 'coolwarm')
        for r in [.5,2.]:ax.add_patch(Circle((0,0),r,fill=False,ls='--',lw=.65,color='black'))
        ax.set(title=title,xlabel='x - center',ylabel='y - center',aspect='equal')
        fig.colorbar(im,ax=ax,label='Surface density' if i<3 else '100 (density / initial density - 1)',extend='both')
    ax = axes[1,2]
    for rows,label,style,width in [(rb,'Bounded','-',2.2),(rp,'Periodic','--',1.1)]:
        x=[.5*(v['lo']+v['hi']) for v in rows]
        ax.semilogy(x,[v['density_abs'] for v in rows],label=label,linestyle=style,linewidth=width)
    ax.axhline(.01,color='black',ls='--',lw=.8,label='1% density reference')
    ax.set(xlabel='Initial radius',ylabel='Mass-weighted mean |density / initial density - 1|',title='Radial density error');ax.legend()
    fig.suptitle('Native W=0 Q1 boundary comparison: common density scale, no smoothing')
    output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=170);plt.close(fig)
    result={'time':mt['t'],'sources':{str(p):digest(p) for p in [initial,endpoint,periodic_initial,periodic_endpoint]},
            'shared_density_limits':[0,1.05],'shared_density_change_percent_limits':[-1,1],
            'identical_initial_disk_IDs_and_physical_data':True,'disk_sites':int(disk.sum()),
            'bounded_radial':rb,'periodic_radial':rp,
            'bounded_disk_density_abs':mean(abs(delta[disk]),a['mass'][disk]),
            'periodic_disk_density_abs':mean(abs(old_delta[original_disk]),c['mass'][original_disk]),
            'scope':'Same-time new bounded versus prior periodic physical evolution with identical scientific disk IC. Exterior/boundary conditions differ. A density image is not a long-time stability proof; colors saturate outside declared fixed limits.'}
    summary.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['initial','endpoint','periodic_initial','periodic_endpoint','output','summary']:
        parser.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    args=parser.parse_args();print(json.dumps(plot(**vars(args)),indent=2))
