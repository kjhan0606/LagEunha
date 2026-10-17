"""Small static pressure-kick diagnostics; no time evolution loop."""
import json
import numpy as np
from test_compatible_pressure_kick import PressureKickTests
from compatible_pressure_kick import kick,fit,face_velocity

maker=PressureKickTests()
records=[]
for dt in (.01,.005):
    states,edges=maker.make(lambda x:.6+.2*x[:,0],lambda x:np.zeros_like(x))
    out,d=kick(states,edges,dt);fits=[fit(s) for s in out]
    jump=max(abs(np.dot(face_velocity(fits[i],(a+b)/2)-face_velocity(fits[j],(a+b)/2),n)) for i,j,a,b,n,L in edges if j>=0)
    records.append(dict(dt=dt,normal_velocity_jump=float(jump),independent_flight_separation_estimate=float(dt*jump),**d))
print(json.dumps(dict(pressure_discontinuity_probe=records),indent=2))
