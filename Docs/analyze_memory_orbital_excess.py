"""Read existing snapshots: excess orbital energy over the circular minimum.

Static canonical central potential only. This is not a conserved hydro energy,
radial action, or proof of physical accuracy. No state or force is modified.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
from audit_canonical_pair_history import gravity_coefficient

def mass(r):
    return gravity_coefficient(r)*r**3

def excess(r, vr, ell, order):
    # All samples here have nonzero angular momentum. No fabricated zero-L limit.
    assert np.all(r>0) and np.all(np.abs(ell)>0)
    lo=np.zeros_like(r); hi=np.maximum(.5,ell**2)
    for _ in range(64):
        mid=(lo+hi)/2
        lower=mass(mid)*mid < ell**2
        lo=np.where(lower,mid,lo); hi=np.where(lower,hi,mid)
    rc=(lo+hi)/2
    root_error=float(np.max(np.abs(mass(rc)*rc-ell**2)/ell**2))
    # Ueff(r)-Ueff(rc) = integral_rc^r [M(s)*s-L^2]/s^3 ds.
    # Avoid subtraction of two O(1) potential energies for nearly circular states.
    nodes,weights=leggauss(order)
    delta=r-rc; integral=np.zeros_like(r)
    for node,weight in zip(nodes,weights):
        s=rc+.5*(node+1)*delta
        integral+=weight*(mass(s)*s-ell**2)/s**3
    result=.5*vr**2+.5*delta*integral
    assert np.isfinite(result).all()
    return result,rc,root_error

# Independent exterior Kepler ellipse identity, different radial phases.
phase=np.linspace(0,2*np.pi,101)
a=2.; ecc=.2
r=a*(1-ecc*np.cos(phase)); ell=np.full_like(r,np.sqrt(a*(1-ecc**2)))
vr=np.sqrt(a)*ecc*np.sin(phase)/r
test,_,rooterr=excess(r,vr,ell,64)
expected=ecc**2/(2*a*(1-ecc**2))
analytic_error=float(np.max(np.abs(test-expected)))
assert analytic_error<1e-13

base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--folder',type=Path,default=base/'Codex_reversible_memory_physical_to_t10_20261018/memory_from_t0')
parser.add_argument('--steps',type=int,nargs='+',default=[562,802,1122])
args=parser.parse_args()
assert all(n>=0 for n in args.steps) and len(set(args.steps))==len(args.steps)
folder=args.folder
env=json.loads((folder/'environment.json').read_text())
assert env['EUNHA_KEPLER_SOFT']=='spline'
for key,value in [('EUNHA_KEPLER_RC',.5),('EUNHA_KEPLER_CX',3),('EUNHA_KEPLER_CY',3)]:
    assert float(env[key])==value
icpath=base/'Codex_canonical_initial_20261017/run/khout.000000.dat'
ic,t0=read(icpath); assert t0==0
r0=np.linalg.norm(np.column_stack((ic['x']-3,ic['y']-3)),axis=1)
masks={'inner_background':(r0<.5)&(ic['den']<.5),
       'inner_shear_band':(r0>=.3)&(r0<.5),
       'disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
rows=[]
for path in [icpath]+[folder/f'khout.{n:06d}.dat' for n in args.steps]:
    data,t=read(path)
    assert np.array_equal(data['id']&np.uint64(2**56-1),ic['id']&np.uint64(2**56-1))
    assert np.array_equal(data['mass'],ic['mass'])
    x=data['x']-3; y=data['y']-3; r=np.hypot(x,y)
    vr=(x*data['vx']+y*data['vy'])/r
    ell=x*data['vy']-y*data['vx']
    e32,_,_=excess(r,vr,ell,32)
    e64,rc,err=excess(r,vr,ell,64)
    row={'snapshot':str(path),'sha256':sha(path),'time':t,
         'circular_radius_relative_root_residual':err,
         'quadrature_max_abs_difference_specific_energy':float(np.max(abs(e32-e64))),
         'minimum_specific_excess_unclipped':float(e64.min()),'regions':{}}
    for name,mask in masks.items():
        row['regions'][name]={'cells':int(mask.sum()),
            'orbital_excess_energy':float(np.sum(data['mass'][mask]*e64[mask])),
            'quadrature_difference_energy':float(np.sum(data['mass'][mask]*abs(e32[mask]-e64[mask]))),
            'radius_minus_circular_radius_mass_RMS':float(np.sqrt(np.average((r[mask]-rc[mask])**2,weights=data['mass'][mask])))}
    rows.append(row)
print(json.dumps({'folder':str(folder),'environment':env,'ellipse_identity_max_absolute_error':analytic_error,'frames':rows,
    'scope':'Eorb-Ecirc(L) removes pure central-gravity epicyclic phase exchange. '
    'Hydro forces and numerical errors can change it; growth alone does not isolate their causes. '
    'Initial material regions, fixed mass, canonical static central potential. No steering.'},indent=2))
