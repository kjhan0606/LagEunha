"""Analytic circular field on saved canonical t0 polygons; diagnostic lift only.
Uniform density, no discontinuous disk interface; not a restart conversion.
"""
import hashlib,json,os,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018')
sys.path.insert(0,'/gpfs/kjhan/LagForce')
import moment_ale_evolve as model
from periodic2d_hllc import quadrature
import quadratic_velocity_reconstruction as quad
from audit_canonical_pair_history import gravity_coefficient

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
rawpath=root/'raw_faces.dat.rhs.rank0000';facepath=root/'raw_faces.dat.rank0000'
raw=np.loadtxt(rawpath);raw=raw[np.argsort(raw[:,0])];x=raw[:,1:3]-3
faces=np.loadtxt(facepath);owner=faces[:,0].astype(int)
active=np.load('/gpfs/kjhan/LagForce/Codex_pressure_rate_backend_v2_20261018/assembly.npz')['active']
ids=np.flatnonzero(active & (np.linalg.norm(x,axis=1)<.7))
def exact(points):
    points=np.atleast_2d(points);r=np.linalg.norm(points,axis=1)
    return np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack([-points[:,1],points[:,0]])
states=[];polys=[];factors=[];mass=[];face_groups=[]
sort=np.argsort(owner,kind='stable');ordered=owner[sort]
for index in ids:
    rows=sort[np.searchsorted(ordered,index):np.searchsorted(ordered,index,side='right')]
    ff=faces[rows];vertices=x[index]+ff[:,6:8]
    angles=np.arctan2(vertices[:,1]-x[index,1],vertices[:,0]-x[index,0]);poly=vertices[np.argsort(angles)]
    pts,dm=quadrature(poly);v=exact(pts);m=dm.sum();q=dm@pts;Q=pts.T@(dm[:,None]*pts)
    p=dm@v;D=pts.T@(dm[:,None]*v);Pi=v.T@(dm[:,None]*v)
    moment=np.r_[m,q,Q.ravel(),p,D.ravel(),Pi.ravel(),m*.002777777777777778/(model.GAMMA-1)+.5*np.trace(Pi)]
    s=model.decode(poly,moment);res=v-s['v']-(pts-s['c'])@s['B'].T
    _,factor=np.linalg.qr(np.sqrt(dm)[:,None]*res,mode='reduced')
    states.append(s);polys.append(poly);factors.append(factor.T);mass.append(m);face_groups.append(ff)
mass=np.array(mass);out=[]
os.environ.update(LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
for mode in ('affine','dominant','projected'):
    for s in states:s.pop('velocity_quadratic',None)
    diag={}
    if mode!='affine':
        os.environ['LAGFORCE_QUADRATIC_BUDGET']=mode
        diag=quad.prepare(states,polys,[],mass,factors)
    sums=np.zeros(2);lengthsum=0.;maxerr=0.;count=0
    gx,gw=np.polynomial.legendre.leggauss(4)
    for local,index in enumerate(ids):
        r=np.linalg.norm(x[index])
        if not .1<r<.45:continue
        for row in face_groups[local]:
            a,b=x[index]+row[6:8],x[index]+row[8:10];length=np.linalg.norm(b-a)
            for z,w in zip(gx,gw):
                point=a+(z+1)/2*(b-a);v=model.velocity(states[local],point);ve=exact(point)[0]
                error=np.linalg.norm(v-ve);sums+=length*w/2*np.array([error**2,ve@ve]);lengthsum+=length*w/2;maxerr=max(maxerr,error)
        count+=1
    out.append(dict(mode=mode,cells=count,velocity_face_RMS=float(np.sqrt(sums[0]/lengthsum)),relative_velocity_RMS=float(np.sqrt(sums[0]/sums[1])),max_velocity_error=maxerr,diagnostic=diag))
print(json.dumps(dict(scope=__doc__,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (rawpath,facepath)},
    reconstruction_cells=len(ids),evaluation_region='.1<generator radius<.45; neighbors from r<.7',results=out,
    limits='Uses analytic circular velocity integrated into new diagnostic moments, not stored production p/D/Pi. Uniform density excludes density-interface errors. One canonical mesh is not a resolution study or physical evolution.'),indent=2))
