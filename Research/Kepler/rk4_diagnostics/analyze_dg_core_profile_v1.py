"""Isolate velocity-curvature and the existing IC's central density cusp.

Offline selected-cell integration only. The original IC is never rewritten.
Constant-density controls are explicit analysis controls, not evolved ICs.
Corner-anchored, geometrically segmented radial quadrature resolves r+1e-6.
"""
from pathlib import Path
import json,struct
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[1]
ids={67943,68000,68527,68584,67935,68001,68526,68592}
gamma=float(np.float32(1.666667));pressure=1.e-6;eps=.047;cx=cy=12.
out=dict(scope='selected-cell independent projection/control analysis; no IC change or evolution',
    regularization_length=1.e-6,epsilon=eps,orders=[],density_samples=[])
for r in (0.,1.e-6,1.e-5,1.e-4,1.e-3,.008286407592029853,.018528970665049098,.047):
    rho=.001+2/(r+1.e-6)*.5*(1+np.tanh((r-2)/.25))*.5*(1+np.tanh((8-r)/.25))
    out['density_samples'].append(dict(r=r,rho=float(rho)))

def load(order):
    case=ROOT/f'rk4_diagnostics/dg_projection_v2_p{order}/fine_fresh'
    path=case/'dg_project_r0.bin'
    with path.open('rb') as f:magic,version,n,p,nb,nm,size,_=struct.unpack('<8Q',f.read(64))
    assert p==order and size==288+40*nb
    dtype=np.dtype(dict(names=['id','x','y','basis','mean','c','metrics'],
        formats=['<u8','<f8','<f8',('<f8',6),('<f8',nb),('<f8',(4,nb)),('<f8',17)],
        offsets=[0,8,16,24,72,72+8*nb,136+40*nb],itemsize=size))
    a=np.fromfile(path,dtype=dtype,offset=64,count=n);rows={int(q['id']):q for q in a if int(q['id']) in ids}
    polygons={}
    with (case/'dg_polygons_r0.bin').open('rb') as f:
        header=struct.unpack('<4Q',f.read(32));assert header[2]==n
        for _ in range(n):
            pid,count,status=struct.unpack('<Qii',f.read(16));data=f.read(16*count)
            if pid in ids:polygons[pid]=np.frombuffer(data,dtype='<f8').reshape(count,2).copy()
    return rows,polygons

def quadrature(poly,basis,anchor,angular_n,radial_n,segmented):
    a,w=leggauss(angular_n);a=(a+1)*.5;w=w*.5
    t,tw=leggauss(radial_n);t=(t+1)*.5;tw=tw*.5
    points=[];weights=[]
    for k in range(len(poly)):
        b,c=poly[k]-anchor,poly[(k+1)%len(poly)]-anchor
        area2=b[0]*c[1]-b[1]*c[0]
        if area2==0:continue
        assert area2>0
        boundary=(1-a[:,None])*b+a[:,None]*c
        breaks=[0.,1.]
        if segmented:
            scale=1.e-6/max(np.linalg.norm(b),np.linalg.norm(c));breaks=[0.,min(scale,1.)]
            while breaks[-1]<1:breaks.append(min(4*breaks[-1],1.))
        for lo,hi in zip(breaks,breaks[1:]):
            tt=lo+(hi-lo)*t
            p=anchor+tt[:,None,None]*boundary[None,:,:]
            ww=area2*np.outer(tw*(hi-lo)*tt,w)
            points.append(p.reshape(-1,2));weights.append(ww.ravel())
    return np.concatenate(points),np.concatenate(weights)

def basis_at(q,row,order):
    b=row['basis'];u=(q[:,0]-b[1])/b[3];v=(q[:,1]-b[2]-b[4]*u)/b[5]
    phi=np.column_stack([u**(d-j)*v**j for d in range(order+1) for j in range(d+1)])
    phi[:,1:]-=row['mean'][1:];return phi

def source(q,row):
    x=row['x']+q[:,0]-cx;y=row['y']+q[:,1]-cy;r=np.hypot(x,y)
    rho=.001+2/(r+1.e-6)*.5*(1+np.tanh((r-2)/.25))*.5*(1+np.tanh((8-r)/.25))
    omega=(r*r+eps*eps)**(-.75);return rho,-y*omega,x*omega

def project(q,w,row,order,check,constant_density=False):
    rho,vx,vy=source(q,row)
    if constant_density:rho[:]=np.sum(w*rho)/row['basis'][0]
    _,v0x,v0y=source(np.zeros((1,2)),row);v0x=float(v0x[0]);v0y=float(v0y[0]);vx-=v0x;vy-=v0y
    u=np.column_stack([rho,rho*vx,rho*vy,pressure/(gamma-1)+.5*rho*(vx*vx+vy*vy)])
    phi=basis_at(q,row,order);integ=np.einsum('i,ij->j',w,u);avg=integ/row['basis'][0]
    h=phi[:,1:].T@(w[:,None]*phi[:,1:]);rhs=phi[:,1:].T@(w[:,None]*(u-avg))
    c=np.zeros((4,phi.shape[1]));c[:,0]=avg;c[:,1:]=np.linalg.solve(h,rhs).T
    v=basis_at(check,row,order)@c.T
    p=(gamma-1)*(v[:,3]-.5*(v[:,1]**2+v[:,2]**2)/v[:,0])
    return p,integ

for order in (2,4,6):
    rows,polygons=load(order);data=dict(order=order,cells=[])
    for pid in sorted(ids):
        row=rows[pid];poly=polygons[pid];anchor=np.array([cx-row['x'],cy-row['y']])
        corner=bool(np.any(np.all(poly==anchor,axis=1)))
        if not corner:anchor=row['basis'][1:3]
        check,cw=quadrature(poly,row['basis'],anchor,24,24,corner)
        t=(leggauss(32)[0]+1)*.5;face=[]
        for k in range(len(poly)):face.extend((1-t[:,None])*poly[k]+t[:,None]*poly[(k+1)%len(poly)])
        check=np.vstack([check,np.array(face),poly,row['metrics'][7:9]])
        prev=None;refs=[]
        for an,rn in [(24,12),(48,20),(64,28)]:
            q,w=quadrature(poly,row['basis'],anchor,an,rn,corner)
            p,integ=project(q,w,row,order,check)
            pc,ic=project(q,w,row,order,check,True)
            record=dict(angular_points=an,radial_points_per_segment=rn,
                min_P=float(p.min()),max_P=float(p.max()),max_relative_P_error=float(np.max(abs(p/pressure-1))),
                constant_density_min_P=float(pc.min()),constant_density_max_relative_P_error=float(np.max(abs(pc/pressure-1))),
                mass=float(integ[0]),constant_density_mass=float(ic[0]))
            if prev is not None:record['relative_P_change']=float(np.max(abs(p-prev))/pressure)
            refs.append(record);prev=p
        data['cells'].append(dict(id=pid,r=float(np.hypot(row['x']-cx,row['y']-cy)),
            corner_anchored_segmented_radial=corner,native_max_P_error=float(row['metrics'][3]),reference=refs))
    out['orders'].append(data)
path=ROOT/'rk4_diagnostics/dg_core_profile_v1.json';path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
