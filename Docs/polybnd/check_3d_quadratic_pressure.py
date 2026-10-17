import argparse,json
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from check_3d_coupled_rhs import initial
from kinetic_moment_rhs_3d import rhs
from moment_state_3d import unpack,decode
from pressure_reconstruction_3d import prepare,value


def run(library):
    geo=Geometry3D(library);rng=np.random.default_rng(736);z=(np.arange(3)+.5)/3;grid=np.meshgrid(z,z,z,indexing='ij');sites=np.column_stack([a.ravel() for a in grid])+rng.uniform(-.01,.01,(27,3));weights=np.zeros(27)
    raw,L=initial(geo,sites,weights,lambda p:np.zeros_like(p));gradient=np.array([.1,-.2,.15]);H=np.array([[.2,.1,-.05],[.1,-.1,.03],[-.05,.03,.1]])
    for r in raw:
        m,q,Q,*_=unpack(r);r[34]=(m+q@gradient+.5*np.sum(Q*H))/(5/3-1)
    dr,gen,z,diag=rhs(geo,sites,weights,raw,L,face_order=3,pressure_reconstruction='quadratic');exact=np.zeros_like(raw)
    for i,r in enumerate(raw):
        m,q,Q,p,D,Pi,E=unpack(r);exact[i,13:16]=-m*gradient-H@q;exact[i,16:25]=(-np.outer(q,gradient)-Q@H.T).ravel()
    error=float(np.max(abs(dr-exact)));assert error<2e-12
    polys,_,_=geo.build(sites,weights);states=[decode(poly,r) for poly,r in zip(polys,raw)];metadata,d=prepare(states,polys);mean_error=field_error=0.
    for s,p,meta in zip(states,polys,metadata):
        x=s['density'].points;w=s['density'].weights;computed=value(meta,x);expected=1+x@gradient+.5*np.einsum('ni,ij,nj->n',x,H,x)
        field_error=max(field_error,float(max(abs(computed-expected))));mean_error=max(mean_error,float(abs(w@computed/w.sum()-s['P'])))
    assert field_error<2e-12 and mean_error<2e-13 and d['pressure_limited_cells']==0
    # Discontinuous positive cell means exercise the Bernstein limiter.
    for i,s in enumerate(states):s['P']=1. if i<13 else .001
    limited,ld=prepare(states,polys);minimum=1.;limited_mean=0.
    for s,p,meta in zip(states,polys,limited):
        x=s['density'].points;w=s['density'].weights;minimum=min(minimum,float(value(meta,x).min()),float(value(meta,p).min()));limited_mean=max(limited_mean,float(abs(w@value(meta,x)/w.sum()-s['P'])))
    assert minimum>0 and limited_mean<2e-13 and ld['pressure_limited_cells']>0
    return dict(passed=True,max_raw_rhs_error=error,field_error=field_error,mean_error=mean_error,limited_mean_error=limited_mean,minimum_sampled_pressure=minimum,limited_cells=ld['pressure_limited_cells'],scope='3D quadratic pressure reproduction, exact Euler moment acceleration, mean preservation and tetrahedral positivity limiter')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');a=p.parse_args();result=run(a.library);Path('quadratic_pressure_3d_report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
