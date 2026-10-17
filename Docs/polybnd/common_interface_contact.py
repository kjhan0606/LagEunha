"""Common physical planar contact across grids; independent static basic test."""
import json,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup,D0
from discontinuous_density_audit import clip
from moment_aligned_density_audit import polygon_moments
from predictive_boundary_density import BoundaryDensity
from predictive_coupling_batch import metrics
from predictive_resolution_audit import trace_error
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def initialize(n,ratio,angle):
    _,old=setup('shear_contact',n,ratio);a=old[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],0.)
    theta=np.radians(angle);normal=np.array([np.cos(theta),np.sin(theta)]);offset=normal@np.array([.5,.5])
    u=np.array([.2,0.]);B=np.array([[0.,1.1],[0.,0.]]);rows=[]
    for v in cells:
        geo=polygon_moments(v);cut=clip(v,normal,offset)
        high=np.zeros(7) if cut.is_empty else polygon_moments(np.asarray(cut.exterior.coords[:-1]))
        raw=geo+(ratio-1)*high;m=raw[0];q=raw[1:3];Q=raw[3:].reshape(2,2);bq=B@q
        p=m*u+bq;D=np.outer(q,u)+Q@B.T
        Pi=m*np.outer(u,u)+np.outer(u,bq)+np.outer(bq,u)+B@Q@B.T
        E=.6*geo[0]/(model.GAMMA-1)+.5*np.trace(Pi)
        rows.append(np.r_[raw,p,D.ravel(),Pi.ravel(),E])
    raw=np.asarray(rows);mass=raw[:,0]
    high=clip(D0,normal,offset);low=clip(D0,-normal,-offset)
    refs=[np.asarray(low.exterior.coords[:-1]),np.asarray(high.exterior.coords[:-1])]
    exact=polygon_moments(D0)+(ratio-1)*polygon_moments(refs[1])
    np.testing.assert_allclose(raw[:,:7].sum(0),exact,rtol=1e-12,atol=1e-13)
    y=np.r_[np.column_stack([a[:,:2],np.zeros(len(mass)),raw[:,1:]]).ravel(),np.zeros(4)]
    return mass,y,refs,[1.,ratio]


def run():
    import sys
    representations=(('uniform3',3,False),('uniform4',4,False),('aligned2',2,True)) if '--refine' in sys.argv else (('uniform2',2,False),)
    results=[]
    for n in (3,4):
        for angle in (0.,23.):
            mass,y,refs,phases=initialize(n,100.,angle);a=y[:-4].reshape(-1,20)
            cells,edges=model.geometry(a[:,:2],a[:,2],0.);L=[np.zeros((2,2)) for _ in mass]
            for mode,tag,level,aligned in [(mode,*rep) for mode in ('neighbor','moment') for rep in representations]:
                started=time.monotonic();provider=BoundaryDensity(mode,vector_geometry=True,core_level=level,align_core=aligned)
                try:
                    dy,d=kinetic.rhs(0.,y,mass,covariance_factors=L,face_density_provider=provider,face_segments=64,face_order=4,exact_bin_transport=True)
                    m=metrics(dy,d,mass);face=trace_error(provider,edges,refs,phases)
                    passed=m['mass_flux']<1e-10 and m['closure_moment_error']<1e-8 and m['density_min']>0 and m['closure_volume_cap_excess']<1e-8 and m['pressure_rate']<1e-8 and max(m['affine_rhs_errors'].values())<1e-8 and m['ledger_rate_error']<1e-10
                    row=dict(n=n,angle=angle,mode=mode,passed=bool(passed),status='validated' if passed else 'gate_failure',metrics=m,face_pointwise_L1_relative=face['relative_L1'],priors=provider.last_priors,seconds=time.monotonic()-started)
                except Exception as exc:row=dict(n=n,angle=angle,mode=mode,passed=False,status='closure_unavailable',error=repr(exc),seconds=time.monotonic()-started)
                row.update(representation=tag,core_level=level,align_core=aligned)
                results.append(row);atomic('common_interface_progress.json',dict(results=results));print(n,angle,mode,tag,row['status'],row.get('face_pointwise_L1_relative'),flush=True)
    # Each of the four common physical initial states needs at least one
    # passing predictive closure. Failed representations remain explicit.
    usable=[dict(mode=mode,representation=tag) for mode in ('neighbor','moment') for tag,_,_ in representations
        if all(any(r['n']==n and r['angle']==angle and r['mode']==mode and r['representation']==tag and r['passed'] for r in results) for n in (3,4) for angle in (0.,23.))]
    report=dict(completed=True,passed=bool(usable),usable_configurations=usable,all_representations_passed=all(r['passed'] for r in results),results=results,
        note='Same exact initial rho=1/100 planar interface across N3/N4. Exact m,q,Q initialization; affine shear velocity and constant pressure. No evolved field injection. References only score reconstruction after RHS. Static test, not spatial evolution convergence or Kepler.')
    atomic('common_interface_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
