"""Identical acoustic initial moments scored with two explicit pressure models."""
import argparse,json,time
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from moment_density_3d import quadrature
from moment_state_3d import pack
from check_3d_acoustic_evolution import reference,score,GAMMA


def run(library,resolutions):
    geo=Geometry3D(library);rows=[];start=time.monotonic()
    for n in resolutions:
        rng=np.random.default_rng(1021+n);grid=np.meshgrid(*([(np.arange(n)+.5)/n]*3),indexing='ij');sites=np.column_stack([a.ravel() for a in grid])+rng.uniform(-.05,.05,(n**3,3))/n;weights=np.zeros(len(sites));polys,_,_=geo.build(sites,weights);raw=[]
        for poly in polys:
            x,w=quadrature(poly,8);rho,P,v=reference(x,0.);dm=w*rho
            raw.append(pack(dm.sum(),dm@x,x.T@(dm[:,None]*x),np.zeros(3),np.zeros((3,3)),np.zeros((3,3)),w@P/(GAMMA-1)))
        raw=np.array(raw)
        for method in ('linear','quadratic'):
            result=score(geo,sites,weights,raw,0.,method);row=dict(n=n,pressure_model=method,score=result,scope='initial reconstruction only, no evolved solution');rows.append(row);print(json.dumps(row),flush=True)
            Path('initial_pressure_3d_audit.json').write_text(json.dumps(dict(rows=rows,seconds=time.monotonic()-start,completed=False),indent=2))
    Path('initial_pressure_3d_audit.json').write_text(json.dumps(dict(rows=rows,seconds=time.monotonic()-start,completed=True),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');p.add_argument('--n',type=int,nargs='+',default=[3,4,6,8]);a=p.parse_args();run(a.library,a.n)
