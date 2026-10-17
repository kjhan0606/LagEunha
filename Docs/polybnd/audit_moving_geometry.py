"""Small irregular-mesh geometry quadrature audit; no fluid evolution."""
import json
import numpy as np
from moving_power_spacetime import integrate
rows=[]
for seed in (731,732,733):
    x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])+np.random.default_rng(seed).uniform(-.04,.04,(4,2))
    v=.3*np.column_stack([.2+1.1*x[:,1],np.zeros(4)])
    w=np.array([.02,0.,0.,0.]);dw=-2*w
    for order,split in [(2,True),(4,True),(8,True),(16,True),(16,False)]:
        r=integrate(x,v,w,dw,order,split)
        if split and order>=8:
            assert max(r[k] for k in ('area_error','first_moment_error','second_moment_error'))<5e-13
        rows.append(dict(seed=seed,**r))
print(json.dumps(rows,indent=2))
