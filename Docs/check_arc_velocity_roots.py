"""Bounded polynomial/root fixtures; execute on compute host."""
import sys,json
from pathlib import Path
import numpy as np
sys.path[:0]=[sys.argv[1],'/gpfs/kjhan/LagForce/Codex_compiled_midpoint_modules_20261018','/gpfs/kjhan/LagForce/Codex_quadratic_thermal_audit_modules_20261018','/gpfs/kjhan/LagForce/Codex_quadratic_split_modules_20261018','/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018','/gpfs/kjhan/LagForce']
from arc_velocity_roots import real_crossings,arc_speed_cuts
from types import SimpleNamespace
from numpy.polynomial import polynomial as P
cases=[[-.5,-.5,.25],[-1.,1.],[-.75,-.1,.2,.85],[-.4,-.4,.3],[-1e-10,1e-10,.7],[0.,0.,0.],[-.99,-.8,-.2,.1,.6,.98]]
rows=[]
for known in cases:
 c=P.polyfromroots(np.asarray(known,dtype=np.longdouble));roots=real_crossings(c)
 for r in roots:assert abs(P.polyval(r,c))<1e-13
 # Even multiplicity roots need not split a non-sign-changing max() kink.
 for r in set(known):
  if known.count(r)%2 and -1<r<1:assert min(abs(r-x) for x in roots)<1e-12
 rows.append(dict(known=known,roots=roots))
state=dict(c=np.array([0.,0.]),v=np.array([1.,.2]),B=np.zeros((2,2)))
arc=SimpleNamespace(radius=2.5,lo=0.,hi=2*np.pi)
cuts=arc_speed_cuts(arc,state,np.zeros((2,2)),1.)
for r in (np.pi/2+np.arctan(.2),3*np.pi/2+np.arctan(.2)):assert min(abs(r-x) for x in cuts)<1e-14
assert real_crossings(np.zeros(7))==[]
print(json.dumps(dict(status='passed_arc_polynomial_fixtures',polynomials=rows,constant_velocity_full_circle_cuts=cuts,provider_calls=0,evolution_steps=0)))
