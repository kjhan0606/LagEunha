"""Independent CPU projection on selected exact native polygons; no evolution.

Uses a Galilean frame in the reference integration to separate polynomial
approximation from cancellation of bulk kinetic energy. This is an equivalent
linear projection of Euler conserved fields, not removal of kinetic energy.
"""
from pathlib import Path
import argparse,json,struct,re
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[1]
pa=argparse.ArgumentParser();pa.add_argument('--order',type=int,default=2);pa.add_argument('--version',default='v2');args=pa.parse_args()
dest=ROOT/f'rk4_diagnostics/dg_projection_{args.version}_p{args.order}'
case=dest/'fine_fresh';path=case/'dg_project_r0.bin'
with path.open('rb') as f:magic,version,n,order,nb,nm,size,rule=struct.unpack('<8Q',f.read(64))
assert magic==0x4746534447503031 and order==args.order and size==288+40*nb
dt=np.dtype(dict(names=['id','x','y','basis','mean','c','integral','coarse_integral','metrics','flags'],
    formats=['<u8','<f8','<f8',('<f8',6),('<f8',nb),('<f8',(4,nb)),('<f8',4),('<f8',4),('<f8',17),('<i4',4)],
    offsets=[0,8,16,24,72,72+8*nb,72+40*nb,104+40*nb,136+40*nb,272+40*nb],itemsize=size))
a=np.fromfile(path,dtype=dt,offset=64,count=n);a.sort(order='id')
assert np.all(a['flags'][:,0]==1)
env=json.loads((case/'environment.json').read_text())
gamma_text=re.search(r'define\s+Gamma_hydro\s*=\s*([0-9.eE+-]+)',(case/'params.dat').read_text()).group(1)
gamma=float(np.float32(float(gamma_text)))
cx=float(env['EUNHA_KEPLER_CX']);cy=float(env['EUNHA_KEPLER_CY']);eps=float(env['EUNHA_KEPLER_EPS'])
rin=float(env['EUNHA_KEPLER_RIN']);rout=float(env['EUNHA_KEPLER_ROUT']);pressure=1.e-6
r=np.hypot(a['x']-cx,a['y']-cy)
indices=set(np.argsort(a['metrics'][:,3])[-8:].tolist())|set(np.argsort(a['metrics'][:,6])[-8:].tolist())
for lo,hi in [(0,.047),(.047,.1),(.1,.2),(.2,.5),(.5,1),(1,2),(2,4),(4,8)]:
    sel=np.flatnonzero((r>=lo)&(r<hi))
    if len(sel):indices.add(int(sel[np.argmax(a['metrics'][sel,3])]))
for old_id in (67935,68001,68526,68592):
    sel=np.flatnonzero(a['id']==old_id)
    if len(sel):indices.add(int(sel[0]))
wanted={int(a['id'][i]) for i in indices};polygons={}
with (case/'dg_polygons_r0.bin').open('rb') as f:
    pmagic,pver,pn,_=struct.unpack('<4Q',f.read(32));assert pmagic==0x4746534447504f31 and pn==n
    for _ in range(pn):
        pid,count,status=struct.unpack('<Qii',f.read(16));data=f.read(16*count)
        if pid in wanted:
            assert status==1;polygons[pid]=np.frombuffer(data,dtype='<f8').reshape(count,2).copy()
assert set(polygons)==wanted

def field(x,y):
    dx=x-cx;dy=y-cy;rr=np.hypot(dx,dy)
    rho=.001+rin/(rr+1.e-6)*.5*(1+np.tanh((rr-rin)/.25))*.5*(1+np.tanh((rout-rr)/.25))
    omega=(rr*rr+eps*eps)**(-.75)
    return rho,-dy*omega,dx*omega

def points(poly,basis,ng):
    node,w=leggauss(ng);node=(node+1)*.5;w=w*.5
    uu,vv=np.meshgrid(node,node,indexing='ij');ww=np.outer(w,w)*(1-uu)
    points=[];weights=[];g=basis[1:3]
    for k in range(len(poly)):
        b,c=poly[k],poly[(k+1)%len(poly)]
        if np.array_equal(b,c):continue
        area2=(b[0]-g[0])*(c[1]-g[1])-(b[1]-g[1])*(c[0]-g[0])
        assert area2>0
        p=(1-uu)[...,None]*(1-vv)[...,None]*g+uu[...,None]*b+(1-uu)[...,None]*vv[...,None]*c
        points.append(p.reshape(-1,2));weights.append((area2*ww).ravel())
    return np.concatenate(points),np.concatenate(weights)

def basis_at(xy,basis,mean):
    u=(xy[:,0]-basis[1])/basis[3];v=(xy[:,1]-basis[2]-basis[4]*u)/basis[5]
    phi=np.column_stack([u**(d-j)*v**j for d in range(order+1) for j in range(d+1)])
    phi[:,1:]-=mean[1:];return phi

def pressure_at(coeff,phi):
    u=phi@coeff.T;return (gamma-1)*(u[:,3]-.5*(u[:,1]**2+u[:,2]**2)/u[:,0])

result=dict(order=order,selected_cells=len(indices),gamma=gamma,
    scope='independent selected-cell quadrature/projection reference; no global or time-evolution certificate',cells=[])
for i in sorted(indices):
    row=a[i];pid=int(row['id']);poly=polygons[pid];basis=row['basis'];mean=row['mean']
    _,v0x,v0y=field(np.array([row['x']]),np.array([row['y']]))
    v0x=float(v0x[0]);v0y=float(v0y[0]);v02=v0x*v0x+v0y*v0y
    # Independent dense assessment grid includes every native face/vertex
    # and its reported worst point, plus a denser interior Gauss rule.
    check,checkw=points(poly,basis,12)
    face=[];t=(leggauss(16)[0]+1)*.5
    for j in range(len(poly)):
        face.extend((1-t[:,None])*poly[j]+t[:,None]*poly[(j+1)%len(poly)])
    check=np.vstack([check,np.array(face),poly,row['metrics'][7:9]])
    phi=basis_at(check,basis,mean)
    native=row['c'].copy()
    native[3]=native[3]-v0x*native[1]-v0y*native[2]+.5*v02*native[0]
    native[1]-=v0x*native[0];native[2]-=v0y*native[0]
    native_p=pressure_at(native,phi)
    refs=[];previous=None
    for ng in (16,32,64):
        q,w=points(poly,basis,ng);ph=basis_at(q,basis,mean)
        rho,vx,vy=field(row['x']+q[:,0],row['y']+q[:,1]);vx-=v0x;vy-=v0y
        u=np.column_stack([rho,rho*vx,rho*vy,pressure/(gamma-1)+.5*rho*(vx*vx+vy*vy)])
        integral=np.einsum('i,ij->j',w,u);average=integral/basis[0]
        h=ph[:,1:].T@(w[:,None]*ph[:,1:])
        rhs=ph[:,1:].T@(w[:,None]*(u-average))
        coef=np.zeros((4,nb));coef[:,0]=average;coef[:,1:]=np.linalg.solve(h,rhs).T
        pp=pressure_at(coef,phi)
        rr=dict(gauss_per_axis=ng,min_P=float(pp.min()),max_P=float(pp.max()),max_P_error=float(np.max(abs(pp/pressure-1))),
            max_native_P_difference=float(np.max(abs(pp-native_p))/pressure),
            mass=float(integral[0]),energy_local=float(integral[3]),gram_condition=float(np.linalg.cond(h)))
        if previous is not None:rr['max_P_change_from_previous']=float(np.max(abs(pp-previous))/pressure)
        previous=pp;refs.append(rr)
    result['cells'].append(dict(id=pid,r=float(r[i]),native_max_P_error=float(row['metrics'][3]),
        native_quad_P_change=float(row['metrics'][6]),reference=refs))
(dest/'reference.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
