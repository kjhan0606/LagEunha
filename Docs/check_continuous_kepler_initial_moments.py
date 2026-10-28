"""Project configured continuous t0 field onto original cells; not a restart.
Measure quadrature convergence including the density discontinuity. No evolution.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np

base = Path('/gpfs/kjhan/LagForce')
sys.path.insert(0, str(base))
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
spec = importlib.util.spec_from_file_location('geom', base/'Codex_weight_remap_moments_centered_20261017/audit.py')
geom = importlib.util.module_from_spec(spec); spec.loader.exec_module(geom)
cache = base/'Codex_pair_hg_modes_t0_20261018'
manifest = json.loads((cache/'manifest.json').read_text())
source = Path(manifest['source'])
assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest['source_sha256']
ic, t = read(source); assert t == 0
env = json.loads((base/'Codex_canonical_initial_20261017/run/environment.json').read_text())
assert env['EUNHA_KEPLER_PROFILE']=='pakmor' and env['EUNHA_KEPLER_BG_ROTATE']=='1'
assert float(env['EUNHA_KEPLER_RC'])==.5
floor=float(env['EUNHA_KEPLER_FLOOR']); rin=float(env['EUNHA_KEPLER_RIN']); rout=float(env['EUNHA_KEPLER_ROUT'])
x=np.column_stack([ic['x']-3,ic['y']-3]); ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7)
v0=np.column_stack([ic['vx'],ic['vy']]); omega=np.sqrt(gravity_coefficient(np.linalg.norm(x[ids],axis=1)))
assert np.max(abs(v0[ids]-omega[:,None]*np.column_stack([-x[ids,1],x[ids,0]])))<1e-12
faces=np.loadtxt(cache/'raw_faces.dat.rank0000'); order=np.argsort(faces[:,0],kind='stable'); owner=faces[order,0]
polys=[]; crossing=[]
for i in ids:
    rows=faces[order[np.searchsorted(owner,i):np.searchsorted(owner,i,side='right')]]
    vertices=rows[:,6:8]; angles=np.arctan2(vertices[:,1],vertices[:,0]); p=vertices[np.argsort(angles)]+x[i]
    polys.append(p)
    # Exact edge minimum radius (these cells do not surround origin at rin).
    a=p; d=np.roll(p,-1,axis=0)-p
    u=np.clip(-np.sum(a*d,axis=1)/np.sum(d*d,axis=1),0,1)
    rmin=np.min(np.linalg.norm(a+u[:,None]*d,axis=1)); rmax=np.max(np.linalg.norm(p,axis=1))
    crossing.append(rmin<rin<rmax)
crossing=np.array(crossing); gamma=float(np.float32(1.666667)); pressure=float(ic['pressure'][0])
assert np.max(abs(ic['pressure']-pressure))<1e-15
root=base/'Codex_continuous_kepler_initial_moments_20261018'; root.mkdir(exist_ok=False)
results=[]; previous=None
for nq in (4,8,16,32):
    moments=[]; internal=[]; stable_thermal=[]; volume_error=[]
    for i,p in zip(ids,polys):
        pts,dv=geom.quadrature(p,nq); r=np.linalg.norm(pts,axis=1)
        rho=np.where((r>=rin)&(r<=rout),1.,floor); dm=rho*dv
        vel=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack([-pts[:,1],pts[:,0]])
        m=dm.sum(); q=dm@pts; Q=pts.T@(dm[:,None]*pts)
        momentum=dm@vel; D=pts.T@(dm[:,None]*vel); Pi=vel.T@(dm[:,None]*vel)
        U=pressure*dv.sum()/(gamma-1); E=U+.5*np.trace(Pi); phi=dm@geom.potential(pts)
        moments.append(np.r_[m,q,Q[0,0],Q[0,1],Q[1,1],momentum,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],E,phi])
        delta=vel-momentum/m; internal.append(.5*np.sum(dm*np.sum(delta*delta,axis=1)))
        stable_thermal.append(U); volume_error.append(dv.sum()/ic['volume'][i]-1)
    q=np.array(moments); assert q.shape==(len(ids),17)
    np.savez(root/f'moments_q{nq}.npz',ids=ids,moments=q,internal_kinetic=internal,thermal=stable_thermal,crosses_inner_edge=crossing)
    summary=dict(order=nq,cells=len(ids),crossing_cells=int(crossing.sum()),minimum_thermal=min(stable_thermal),
        max_volume_relative_error=float(np.max(abs(np.array(volume_error)))),
        mass_difference_from_point=float(np.sum(q[:,0]-ic['mass'][ids])),
        max_cell_relative_mass_difference=float(np.max(abs(q[:,0]/ic['mass'][ids]-1))),
        angular_difference_from_point=float(np.sum(q[:,9]-q[:,10]-ic['mass'][ids]*(x[ids,0]*v0[ids,1]-x[ids,1]*v0[ids,0]))),
        total_energy_difference_from_point=float(np.sum(q[:,15]+q[:,16]-(ic['ie'][ids]+.5*ic['mass'][ids]*np.sum(v0[ids]**2,axis=1)+ic['mass'][ids]*geom.potential(x[ids])))),
        sum_internal_kinetic=float(np.sum(internal)))
    if previous is not None:
        summary['successive_relative_mass_max']={name:float(np.max(abs(q[mask,0]/previous[mask,0]-1))) for name,mask in [('crossing',crossing),('noncrossing',~crossing)]}
        summary['successive_moment_max_abs']=float(np.max(abs(q-previous)))
    previous=q;results.append(summary)
report=dict(scope=__doc__,source=manifest,environment=env,results=results,
    limits='Finite patch r(generator)<.7. Continuous radial density differs from point-sampled IC mass. Positive thermal energy is constructed, not earned by evolution. Unsplitted circle interface quadrature must converge before use. No force, restart, CUDA/MPI, or stability pass.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
