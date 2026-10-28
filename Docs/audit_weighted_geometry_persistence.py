"""Prescribed orbit replay of saved static weighted-geometry candidates."""
import argparse
import importlib.util
import json
from pathlib import Path
import time

import numpy as np


def main():
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);args=ap.parse_args()
    args.output.mkdir(exist_ok=False)
    prior=Path('/gpfs/kjhan/LagForce/Codex_weighted_compression_geometry_20261017')
    spec=importlib.util.spec_from_file_location('weighted_geometry',prior/'audit.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator')
    a=np.loadtxt(source/'rhs.dat');a=a[np.argsort(a[:,0])]
    x=a[:,1:3]-3;v=a[:,3:5];r=np.linalg.norm(x,axis=1)
    radial=np.sum(x*v,axis=1)/r
    assert np.max(abs(radial[r<.7]))<1e-12
    omega=(x[:,0]*v[:,1]-x[:,1]*v[:,0])/r**2
    omega[r>=2.5]=0  # outer rotation is irrelevant to measured central patch
    ids=np.flatnonzero(r<.7);check=set(ids[::max(1,len(ids)//16)].tolist())
    data=np.load(prior/'damp1_fraction1.npz');w=data['weights']
    weight_volume_scale=float(np.max(abs(w)/a[:,6]))
    records=[];started=time.monotonic()
    for t in (0.,.01,.05,.1,.25,.5,1.):
        c,s=np.cos(omega*t),np.sin(omega*t)
        pos=np.column_stack((c*x[:,0]-s*x[:,1],s*x[:,0]+c*x[:,1]))
        vel=np.column_stack((-omega*pos[:,1],omega*pos[:,0]))
        for name,weights in (('voronoi',np.zeros(len(x))),('constant_saved_weights',w)):
            try:
                vol,rate,L,cert=module.geometry(pos,weights,vel,ids,check)
                inner=r[ids]<.5;guard=~inner
                rec={'elapsed':t,'candidate':name,'valid_geometry':True,'sample_plane_violation':cert,
                     'inner_generator_volume_rate_rms':float(np.sqrt(np.average((rate[inner]/vol[inner])**2,weights=a[ids[inner],6]))),
                     'guard_generator_volume_rate_rms':float(np.sqrt(np.average((rate[guard]/vol[guard])**2,weights=a[ids[guard],6]))),
                     'min_volume_over_source':float(np.min(vol/a[ids,6])),'max_volume_over_source':float(np.max(vol/a[ids,6]))}
            except (RuntimeError,AssertionError) as exc:
                rec={'elapsed':t,'candidate':name,'valid_geometry':False,'error':repr(exc)}
            records.append(rec);print(json.dumps(rec),flush=True)
    out={'scope':'Kinematic prescribed-circle replay from actual t2 circular-state positions. Saved weights remain constant; generator trajectories identical between candidates. No fluid, mass, pressure or gravity evolution is simulated.',
         'initial_inner_radial_velocity_max_abs':float(np.max(abs(radial[r<.7]))),
         'max_weight_over_source_volume':weight_volume_scale,'records':records,'seconds':time.monotonic()-started,
         'interpretation':'A static improvement that disappears during replay needs an evolving weight law. This replay alone cannot verify conservative mass coupling or long-time Kepler stability.'}
    (args.output/'result.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
