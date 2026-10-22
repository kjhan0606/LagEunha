"""Integrate the same continuous IC over true periodic Voronoi cells.

Only initial data and quadrature are computed here. Evolution remains native.
"""
from pathlib import Path
import hashlib,json,struct,sys,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.spatial import Voronoi
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
dest=root/'core_refinement';dest.mkdir(exist_ok=False)
layout,_=read_snapshot(root/'linear_velocity_limiter/physics_v4/sedov_vector_raw/khout.000001.dat')
h=24/256;eps=.047;gamma=float(np.float32(1.666667));P=float(np.float32(1e-6))
ids=np.arange(256**2);coarse=np.column_stack(((ids%256+.5)*h,(ids//256+.5)*h))
selected=np.hypot(coarse[:,0]-12,coarse[:,1]-12)<.5
assert selected.sum()==88
sub=np.array([(a,b) for b in range(8) for a in range(8)])
refined=coarse[selected,None,:]+((sub[None,:,:]+.5)/8-.5)*h
fine=np.vstack((coarse[~selected],refined.reshape(-1,2)))
parents=np.r_[ids[~selected],np.repeat(ids[selected],64)]
logf0=-np.logaddexp(0.,16.)-np.logaddexp(0.,-64.)
f0=np.exp(logf0)
def shape(r):
    return np.exp(-np.logaddexp(0.,8*(2-r))-np.logaddexp(0.,8*(r-8)))
def density(r):return .001+2*shape(r)/(r+1e-6)

def geometry(xy):
    # All 3x3 periodic images, rather than an unproved boundary strip.
    blocks=[xy]+[xy+np.array([dx,dy]) for dx in (-24,0,24) for dy in (-24,0,24) if dx or dy]
    vor=Voronoi(np.vstack(blocks));volume=np.empty(len(xy));radius=np.empty(len(xy))
    apex=[];ends_a=[];ends_b=[];cell=[];special=[]
    for i in range(len(xy)):
        region=vor.regions[vor.point_region[i]]
        if not region or -1 in region:raise ValueError('unbounded periodic cell')
        v=vor.vertices[region]-12
        cross=v[:,0]*np.roll(v[:,1],-1)-v[:,1]*np.roll(v[:,0],-1)
        local=v-(xy[i]-12)
        area2=np.sum(local[:,0]*np.roll(local[:,1],-1)-local[:,1]*np.roll(local[:,0],-1))
        if area2<0:v=v[::-1];local=local[::-1];area2=-area2
        volume[i]=.5*area2;radius[i]=np.max(np.linalg.norm(local,axis=1))
        # The closest omitted periodic images are at least24 away.
        if not 2*radius[i]<24:raise ValueError('periodic completeness proof failed')
        nxt=np.roll(v,-1,axis=0);edge=nxt-v
        signed=edge[:,0]*(-v[:,1])-edge[:,1]*(-v[:,0])
        roundoff=128*np.finfo(float).eps*np.maximum(1,np.max(abs(v))**2)
        origin_inside=bool(np.all(signed>=-roundoff))
        origin=np.zeros(2) if origin_inside else xy[i]-12
        for a,b in zip(v,nxt):
            c=np.cross(a-origin,b-origin).item()
            if c< -roundoff:raise ValueError('nonconvex triangle fan')
            if c<=0:continue
            apex.append(origin);ends_a.append(a);ends_b.append(b);cell.append(i);special.append(origin_inside)
    if not np.all(volume>0) or abs(volume.sum()-576)>1e-8:raise ValueError('invalid periodic area sum')
    return volume,radius,np.array(apex),np.array(ends_a),np.array(ends_b),np.array(cell),np.array(special)

def integrate(n,geom,selected_cells=None):
    V,R,O,A,B,owner,special=geom
    active=np.ones(len(owner),bool) if selected_cells is None else selected_cells[owner]
    idx=np.flatnonzero(active);result=np.zeros(len(V))
    nodes,weights=leggauss(n);nodes=(nodes+1)/2;weights=weights/2
    u=nodes[None,:,None];v=nodes[None,None,:]
    wu=weights[None,:,None];wv=weights[None,None,:]
    for start in range(0,len(idx),512):
        ix=idx[start:start+512];o=O[ix];a=A[ix]-o;b=B[ix]-o
        edge=(1-v[...,None])*a[:,None,None,:]+v[...,None]*b[:,None,None,:]
        points=o[:,None,None,:]+u[...,None]*edge
        r=np.linalg.norm(points,axis=-1);rho=density(r)
        central=special[ix]
        if central.any():
            rc=r[central]
            logf=-np.logaddexp(0.,8*(2-rc))-np.logaddexp(0.,8*(rc-8))
            rho[central]=.001+2*f0*np.expm1(logf-logf0)/(rc+1e-6)
        q=np.sum(wu*wv*u*rho,axis=(1,2))
        if central.any():
            rr=np.linalg.norm(edge[central,0,:,:],axis=-1)
            # Exact radial integral of the mild1/(r+1e-6) central term.
            analytic=2*f0*(rr-1e-6*np.log1p(rr/1e-6))/(rr*rr)
            q[central]+=np.sum(weights[None,:]*analytic,axis=1)
        cross=a[:,0]*b[:,1]-a[:,1]*b[:,0]
        result+=np.bincount(owner[ix],weights=cross*q,minlength=len(V))
    return result

manifest={'purpose':'fixed-eps initial core resolution discrimination; not long-term stability',
          'eps':eps,'coarse_spacing':h,'fine_spacing':h/8,'coarse_cells_split':int(selected.sum()),
          'continuous_density':'.001+2/(r+1e-6)*[1+tanh((r-2)/.25)]/2*[1+tanh((8-r)/.25)]/2',
          'quadrature':'actual periodic Voronoi polygons; converged triangle Gauss quadrature, analytic central singular radial term',
          'velocity':'original pointwise Plummer circular velocity at generators; no analytic velocity forcing during evolution',
          'output_policy':'two ICs about79MiB total; later native checkpoints only as scientifically required, no repeated whole reruns',
          'cases':{}}
for label,xy,parent in (('coarse_integrated',coarse,ids),('fine_integrated',fine,parents)):
    start=time.monotonic();geom=geometry(xy);V=geom[0]
    m4=integrate(4,geom);mass=integrate(8,geom)
    err=abs(mass-m4)/np.maximum(mass,1e-300)
    refinements=[]
    for order in (16,32,64):
        bad=err>1e-8
        if not bad.any():break
        revised=integrate(order,geom,bad)
        err[bad]=abs(revised[bad]-mass[bad])/np.maximum(revised[bad],1e-300)
        mass[bad]=revised[bad]
        refinements.append({'order':order,'cells':int(bad.sum())})
    if np.max(err)>1e-8:raise ValueError(f'{label}: quadrature has not converged: {err.max()}')
    raw=np.zeros(len(xy)*layout.dtype.itemsize,dtype=np.uint8)
    p=np.ndarray(len(xy),dtype=layout.dtype,buffer=raw)
    r=np.linalg.norm(xy-12,axis=1);omega=(r*r+eps*eps)**(-.75)
    p['u4if']=np.arange(len(xy),dtype=np.uint64);p['x']=xy[:,0];p['y']=xy[:,1]
    p['mass']=mass;p['volume']=V;p['den']=mass/V;p['pressure']=P;p['ie']=P*V/(gamma-1)
    p['vx']=-(xy[:,1]-12)*omega;p['vy']=(xy[:,0]-12)*omega
    d=dest/label;d.mkdir();path=d/'khout.000001.dat'
    with path.open('wb') as f:
        f.write(struct.pack('<4i',-1,len(p),layout.dtype.itemsize,1));f.write(raw.tobytes())
        f.write(struct.pack('<ddfiiddff',0.,1e-7,1e-7,256,256,24.,24.,0.,0.))
    params=(root/'linear_velocity_limiter/cold_matrix/legacy_raw/params.dat').read_text()
    (d/'params.dat').write_text(params)
    np.savez_compressed(d/'geometry.npz',xy=xy,parent=parent,volume=V,quad_error=err)
    row={'n':len(p),'mass':float(mass.sum()),'internal_energy':float(p['ie'].astype(float).sum()),
         'volume_sum':float(V.sum()),'min_volume':float(V.min()),'minimum_generator_radius':float(r.min()),
         'max_mass_quadrature_relative_change':float(err.max()),'quadrature_refinements':refinements,
         'density_cell_average_min':float((mass/V).min()),'density_cell_average_max':float((mass/V).max()),
         'preparation_seconds':time.monotonic()-start,'IC_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    manifest['cases'][label]=row
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(label,row,flush=True)
a,b=manifest['cases'].values()
rel=abs(a['mass']-b['mass'])/a['mass']
if rel>1e-8:raise ValueError(f'continuous total mass mismatch {rel}')
manifest['coarse_fine_total_mass_relative_difference']=rel
(dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
