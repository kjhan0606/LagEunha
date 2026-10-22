"""Read native KH snapshots using compiler-reported offsets; audit exact W=0 cell areas."""
from pathlib import Path
import argparse, json, struct, time
import numpy as np
from scipy.spatial import Voronoi

ROOT=Path(__file__).resolve().parent

def read_snapshot(path):
    path=Path(path)
    with path.open('rb') as f:
        magic,n,psize,stress=struct.unpack('<4i',f.read(16))
        if magic != -1: raise ValueError('Only the native full-stride format is supported')
        f.seek(16+n*psize)
        meta=struct.unpack('<ddfiiddff',f.read(52))
    layout=json.loads((ROOT/'snapshot_offsets.json').read_text())
    if psize != layout['particle_size']: raise ValueError('Compiled layout differs from snapshot stride')
    fields=layout['fields']
    dtype=np.dtype({'names':list(fields),'formats':['<u8' if k=='u4if' else '<f'+str(v['size']) for k,v in fields.items()], 'offsets':[v['offset'] for v in fields.values()],'itemsize':psize})
    data=np.memmap(path,mode='r',dtype=dtype,offset=16,shape=(n,))
    return data,dict(t=meta[0],dt=meta[1],dtold=meta[2],nx=meta[3],ny=meta[4],Lx=meta[5],Ly=meta[6],alpha_vis=meta[7],beta_vis=meta[8],n=n,psize=psize)

def audit(path,cx=12.,cy=12.,eps=.047):
    p,m=read_snapshot(path)
    if np.any(p['w2'] != 0): raise ValueError('This exact area audit is W=0 only')
    xy=np.column_stack([p['x'],p['y']])
    vor=Voronoi(xy)
    area=np.full(len(p),np.nan)
    radius=np.full(len(p),np.nan)
    for i,reg in enumerate(vor.point_region):
        ind=vor.regions[reg]
        if not ind or -1 in ind: continue
        v=vor.vertices[ind]-xy[i]
        area[i]=.5*abs(np.dot(v[:,0],np.roll(v[:,1],-1))-np.dot(v[:,1],np.roll(v[:,0],-1)))
        radius[i]=np.sqrt(np.max(np.sum(v*v,axis=1)))
    edge=np.minimum.reduce([p['x'],p['y'],m['Lx']-p['x'],m['Ly']-p['y']])
    # A 2*circumradius neighborhood inside the original box excludes every
    # omitted periodic-image generator, so these are complete periodic cells.
    safe=(edge>2*radius)&np.isfinite(area)&(area>0)
    dx=p['x']-cx; dy=p['y']-cy; r=np.hypot(dx,dy)
    vr=(dx*p['vx']+dy*p['vy'])/np.maximum(r,1e-30)
    vphi=(dx*p['vy']-dy*p['vx'])/np.maximum(r,1e-30)
    vc=r*(r*r+eps*eps)**(-.75)
    rel=abs(p['volume']-area)/area
    density_rel=abs(p['den']-p['mass']/p['volume'])/np.maximum(abs(p['den']),1e-30)
    out={'snapshot':str(Path(path).resolve()),'meta':m,'area_audit':'global unweighted Voronoi; only cells proved independent of missing periodic images', 'safe_cells':int(safe.sum()),'mass':float(p['mass'].sum()),'Lz':float(np.sum(p['mass']*(dx*p['vy']-dy*p['vx']))),'Eint':float(p['ie'].sum()),'Etotal':float(np.sum(p['ie']+.5*p['mass']*(p['vx']**2+p['vy']**2)-p['mass']/np.sqrt(r*r+eps*eps))), 'max_stored_density_mass_volume_relative':float(density_rel.max()),'radial':[]}
    for lo,hi in [(0,.5),(.5,1),(1,2),(2,4),(4,8),(8,10)]:
        sel=(r>=lo)&(r<hi)&safe
        if not sel.any(): continue
        rr=rel[sel]
        out['radial'].append({'Rlo':lo,'Rhi':hi,'n':int(sel.sum()), 'area_rel_max':float(rr.max()), 'area_rel_p99':float(np.quantile(rr,.99)), 'area_rel_median':float(np.median(rr)), 'cells_area_error_gt_1pct':int((rr>.01).sum()), 'L1_vphi':float(np.mean(abs(vphi[sel]-vc[sel])/vc[sel])), 'rms_vr_over_vc':float(np.sqrt(np.mean((vr[sel]/vc[sel])**2))), 'pressure_median':float(np.median(p['pressure'][sel])), 'pressure_max':float(p['pressure'][sel].max())})
    return out

if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('snapshot'); a.add_argument('output'); args=a.parse_args()
    start=time.monotonic(); result=audit(args.snapshot); result['analysis_seconds']=time.monotonic()-start
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
