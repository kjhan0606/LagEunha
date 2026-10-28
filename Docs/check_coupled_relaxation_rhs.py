"""Experimental full moment hydro RHS with opt-in affine contact relaxation.

Instantaneous small shear patch only. Existing donor residual closure is still
an unvalidated research assumption; this is not a complete Euler/Kepler test.
"""
import json
import numpy as np
from adaptive_moment_batch import setup
from compatible_moment_rhs import rhs

results=[]
for amplitude in [0.,.01]:
    mass,y=setup('shear',4,0.);a=y[:-4].reshape(len(mass),20)
    # Per-cell velocity shift with all raw kinetic moments and energy updated.
    for i,m in enumerate(mass):
        dv=amplitude*np.array([(-1.)**i,(-1.)**(i//4)])
        q=a[i,3:5];p=a[i,9:11].copy();D=a[i,11:15].reshape(2,2)
        deltaPi=np.outer(dv,p)+np.outer(p,dv)+m*np.outer(dv,dv)
        a[i,9:11]+=m*dv;D+=np.outer(q,dv)
        a[i,15:19]+=deltaPi.ravel();a[i,19]+=.5*np.trace(deltaPi)
    off,doff=rhs(0.,y,mass,linear_pressure=True)
    explicit,_=rhs(0.,y,mass,linear_pressure=True,relaxation_alpha=0.)
    assert np.array_equal(off,explicit)
    on,don=rhs(0.,y,mass,linear_pressure=True,relaxation_alpha=1.)
    difference=on-off;d=difference[:-4].reshape(len(mass),20)
    assert np.array_equal(d[:,:9],np.zeros_like(d[:,:9]))
    assert np.array_equal(difference[-4:],np.zeros(4))
    thermal=d[:,19]-.5*(d[:,15]+d[:,18])
    angular=np.sum(d[:,12]-d[:,13]);momentum=np.sum(d[:,9:11],axis=0);energy=np.sum(d[:,19])
    Rrate=[]
    for i,m in enumerate(mass):
        c=a[i,3:5]/m;v=a[i,9:11]/m
        S=a[i,5:9].reshape(2,2)-m*np.outer(c,c)
        C=a[i,11:15].reshape(2,2)-m*np.outer(c,v);G=np.linalg.solve(S,C).T
        dc=d[i,11:15].reshape(2,2)-np.outer(c,d[i,9:11])
        dg=np.linalg.solve(S,dc).T
        dT=d[i,15:19].reshape(2,2)-np.outer(v,d[i,9:11])-np.outer(d[i,9:11],v)
        Rrate.append(np.max(abs(dT-dg@S@G.T-G@S@dg.T)))
    row=dict(amplitude=amplitude,heat=don['relaxation_heat_rate'],
             thermal_rate_sum=float(sum(thermal)),minimum_thermal_rate=float(thermal.min()),
             energy_residual=float(energy),angular_residual=float(angular),momentum_residual=float(abs(momentum).max()),
             residual_covariance_rate_max=float(max(Rrate)),
             max_RHS_change=float(abs(difference).max()),default_equals_explicit_off=True,
             contacts=don['relaxation_contacts'])
    assert abs(energy)<1e-12 and abs(angular)<1e-12 and abs(momentum).max()<1e-12
    assert max(Rrate)<1e-12 and thermal.min()>-1e-12
    assert abs(sum(thermal)-row['heat'])<1e-12
    if amplitude==0:assert row['heat']<1e-20
    else:assert row['heat']>0
    results.append(row)
print(json.dumps(dict(scope=__doc__,cases=results,
    limits=['Single RHS, not a time evolution or source of long-term stability evidence.',
    'No boundary damping; no production C/CUDA/MPI caller connected.',
    'alpha=1 is a diagnostic value, not a validated Kepler viscosity.']),indent=2))
