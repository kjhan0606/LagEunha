"""Coarse canonical-physics control of the full moment/gravity/frozen coupling.

The nominal ring spacing is L/n. Small layouts retain the actual unpadded
count and are labelled auxiliary controls; the official192² IC is unchanged.
"""
import argparse,json,time,os
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from kepler_disk_reference import ring_sites,physical_raw,frozen,CENTER,PRESSURE
from kepler_spline_gravity import acceleration,potential
from factor_moment_evolve import rebuild
from gravity_coupled_2d import step,potential_energy


def install():
    base.domain=lambda t:np.array([[0.,0.],[6.,0.],[6.,6.],[0.,6.]])
    base.boundary_velocity=lambda p:np.zeros(2)


def run(n,dt):
    if Path('initial_state.json').exists():raise RuntimeError('existing IC must be reused explicitly, not overwritten')
    install();sites=ring_sites(n,require_count=(n==192));weights=np.zeros(len(sites));polys,_=base.geometry(sites,weights,0.)
    mass=[];rows=[];factors=[];reference_error=0.;volume_error=0.;start=time.monotonic()
    for i,(site,poly) in enumerate(zip(sites,polys)):
        try:low,_,_=physical_raw(poly,16);raw,L,d=physical_raw(poly,24)
        except Exception:
            Path('initial_quadrature_failure.json').write_text(json.dumps(dict(cell=i,polygon=poly.tolist())));raise
        error=float(np.max(abs(raw-low))/max(np.max(abs(raw)),np.finfo(float).tiny));reference_error=max(reference_error,error)
        if error>1e-10:raise RuntimeError('unresolved Kepler IC quadrature at cell '+str(i)+': '+str(error))
        mass.append(raw[0]);rows.append(np.r_[site,0.,raw[1:]]);factors.append(L)
        if i%32==0:print('Kepler IC',i,'/',len(sites),flush=True)
    mass=np.array(mass);original=np.r_[np.array(rows).ravel(),np.zeros(4)];y=rebuild(original,mass,factors);a=y[:-4].reshape(-1,20)
    initial=dict(t=0.,mass=mass.tolist(),state=y.tolist(),factors=[L.tolist() for L in factors],nominal_spacing_n=n,actual_cells=len(mass),official_resolution=n==192,pressure=PRESSURE,reference_quadrature_relative=reference_error,Pi_initial_completion=float(np.max(abs(y-original))))
    Path('initial_state.json').write_text(json.dumps(initial))
    pressures=[];density_fit=0.
    for i,(poly,m,row) in enumerate(zip(polys,mass,a)):
        try:s=base.decode(poly,np.r_[m,row[3:]])
        except Exception:
            Path('initial_decode_failure.json').write_text(json.dumps(dict(cell=i,polygon=poly.tolist(),raw=np.r_[m,row[3:]].tolist())));raise
        pressures.append(s['P']);density_fit=max(density_fit,s['density'].residual)
    pressure_error=float(np.max(abs(np.array(pressures)/PRESSURE-1)));assert pressure_error<1e-8
    phi=lambda p:potential(p,center=CENTER,rc=.5);gravity=lambda p:acceleration(p,center=CENTER,rc=.5)
    initial_total=float(a[:,19].sum())+potential_energy(0.,y,mass,phi);before=y.copy();oldL=[L.copy() for L in factors];mask=frozen(0.,sites)
    Path('initial_audit.json').write_text(json.dumps(dict(reference_quadrature_relative=reference_error,pressure_relative_error=pressure_error,density_fit=density_fit,total_energy=initial_total,frozen_count=int(mask.sum())),indent=2))
    out,nextL,d=step(0.,y,mass,factors,dt,gravity,frozen=frozen)
    state=out[:-4].reshape(-1,20);np.testing.assert_array_equal(state[mask],a[mask])
    for i in np.flatnonzero(mask):np.testing.assert_array_equal(nextL[i],oldL[i])
    final_total=float(state[:,19].sum())+potential_energy(dt,out,mass,phi)
    checkpoint=dict(t=dt,mass=mass.tolist(),state=out.tolist(),factors=[L.tolist() for L in nextL],nominal_spacing_n=n,actual_cells=len(mass),dt=dt)
    Path('checkpoint.json').write_text(json.dumps(checkpoint))
    result=dict(completed=True,scope='first coupled fluid step, not100-orbit validation',official_resolution=n==192,nominal_spacing_n=n,actual_cells=len(mass),dt=dt,frozen_count=int(mask.sum()),frozen_state_exact=True,initial_pressure_relative_error=pressure_error,reference_quadrature_relative=reference_error,total_energy_change=final_total-initial_total,gas_boundary_energy_ledger=float(out[-2]-before[-2]),gravity_diagnostics=d,seconds=time.monotonic()-start)
    Path('kepler_preflight_report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=16);p.add_argument('--dt',type=float,default=1e-5);a=p.parse_args();run(a.n,a.dt)
