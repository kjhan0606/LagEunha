"""Pressure work split on actual t6 geometry; no modified evolution."""
import argparse, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_t6_entropy_force_probes_20261018'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,default=root)
root=parser.parse_args().root
manifest=json.loads((root/'manifest.json').read_text())
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat')
seed,_=read(b/'Codex_coupled_onset_seed_20261018/reference/khout.000184.dat')
g=float(np.float32(1.666667));r0=np.linalg.norm(np.column_stack((ic['x']-3,ic['y']-3)),axis=1)
masks={'inner_background':(r0<.5)&(ic['den']<.5),'inner_shear_band':(r0>=.3)&(r0<.5)}
out={'scope':'Frozen pressure-force/work decomposition. Not causal fractions or complete perturbation energy budget; TCF/HG and gravity/reference work must still be included.', 'cases':{}}
for name in ('adaptive','legacy'):
    d=root/name;assert (d/'exit.code').read_text().strip()=='0'
    a,t=read(d/'khout.001783.dat')
    assert sha(d/'khout.001783.dat')==manifest['cases'][name]['source_sha256']
    s=np.loadtxt(d/'raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
    f=np.loadtxt(d/'raw_faces.dat.common.rank0000');q=np.loadtxt(d/'raw_faces.dat.rhs.faces.rank0000')
    assert np.array_equal(s[:,0],np.arange(len(a)))
    for col,key in ((1,'x'),(2,'y'),(3,'vx'),(4,'vy'),(5,'pressure'),(6,'volume'),(7,'mass')):
        assert np.allclose(s[:,col],a[key],rtol=1e-12,atol=1e-14), (name,key)
    i,j=f[:,:2].astype(int).T;dx=f[:,2:4];A=f[:,4];distance=np.linalg.norm(dx,axis=1);e=dx/distance[:,None]
    tangent=f[:,5:7]-.5*dx;tangent-=np.sum(tangent*e,axis=1)[:,None]*e
    def force(p):
        z=-A[:,None]*(.5*(p[i]+p[j])[:,None]*e+((p[j]-p[i])/distance)[:,None]*tangent)
        return np.column_stack([np.bincount(i,weights=z[:,k],minlength=len(a)) for k in range(2)])
    p0=ic['pressure']*(ic['volume']/s[:,6])**g
    p2=seed['pressure']*(seed['volume']/s[:,6])**g
    parts={k:force(p) for k,p in [('actual',s[:,5]),('initial_entropy',p0),('inherited_entropy',p2-p0),('post_t2_entropy',s[:,5]-p2)]}
    closure=float(np.max(abs(parts['actual']-parts['initial_entropy']-parts['inherited_entropy']-parts['post_t2_entropy'])))
    assert closure<1e-12
    fq=np.zeros((len(a),2));np.add.at(fq,q[:,0].astype(int),-(q[:,4]*q[:,14])[:,None]*q[:,6:8])
    x=s[:,1:3]-3;r=np.linalg.norm(x,axis=1)
    active=r<1.9
    reproduction=float(np.linalg.norm((parts['actual']-s[:,8:10]+fq)[active])/np.linalg.norm((s[:,8:10]-fq)[active]))
    assert reproduction<1e-4
    vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]));dv=s[:,3:5]-vc
    row={'time':t,'source_sha256':sha(d/'khout.001783.dat'),'force_closure_max_abs':closure,'force_reproduction_relative_L2':reproduction,'regions':{}}
    for label,mask in masks.items():
        row['regions'][label]={k:float(np.sum(dv[mask]*v[mask])) for k,v in parts.items()}
    out['cases'][name]=row
path=root/'entropy_pressure_work.json'
if path.exists(): assert json.loads(path.read_text())==out
else:
    with path.open('x') as stream:json.dump(out,stream,indent=2)
print(json.dumps(out,indent=2))
