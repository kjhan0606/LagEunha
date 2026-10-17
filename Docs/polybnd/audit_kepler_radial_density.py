"""Read-only radial localization of density error; no acceptance exclusions.

Fixed bands isolate both density jumps. Their union is the whole r<2.5
scoring domain. This diagnostic neither evolves nor repairs saved states.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from audit_kepler_energy_step import configure,phi,atomic
from kepler_disk_reference import polar_rule,fields,CENTER,frozen
from potential_moment_coupled_2d import provider

EDGES=np.array([0.,.45,.55,1.95,2.05,2.5])


def run(source,configuration,order):
    data=Path(source).read_bytes();cp=json.loads(data)
    configure(json.loads(Path(configuration).read_text()))
    mass=np.asarray(cp['mass']);a=np.asarray(cp['state'])[:-4].reshape(-1,20)
    polys,_=base.geometry(a[:,:2],a[:,2],cp['t'])
    density_provider=None
    if 'potential_moments' in cp:
        density_provider=provider(np.asarray(cp['potential_moments']),phi,
            lambda p:polar_rule(p,16,radii=(.25,.5)),frozen(cp['t'],a[:,:2]))
    total=np.zeros((5,4));minimum=np.full(5,np.inf);maximum=np.zeros(5)
    worst=[None]*5
    for i,(poly,row,m) in enumerate(zip(polys,a,mass)):
        raw=np.r_[m,row[3:]]
        den=base.Density(poly,m,raw[1:3],raw[3:7].reshape(2,2)) if density_provider is None else density_provider(i,poly,raw)
        x,w=polar_rule(poly,order,radii=(.25,.45,.5,.55,1.95,2.,2.05,2.5))
        r=np.linalg.norm(x-CENTER,axis=1);keep=r<2.5;x=x[keep];w=w[keep];r=r[keep]
        if not len(w):continue
        rho=fields(x)[0];rn=den.value(x);bins=np.searchsorted(EDGES,r,side='right')-1
        for j in range(5):
            use=bins==j
            if not np.any(use):continue
            q=x[use];z=rn[use];ref=rho[use];ww=w[use]
            total[j]+=ww@np.column_stack([np.ones(len(ww)),ref,z,abs(z-ref)])
            minimum[j]=min(minimum[j],float(z.min()))
            rel=abs(z/ref-1);k=int(np.argmax(rel))
            if worst[j] is None or rel[k]>maximum[j]:
                maximum[j]=float(rel[k]);worst[j]=dict(cell=i,point=q[k].tolist(),density=float(z[k]),reference=float(ref[k]))
    rows=[dict(r_min=float(EDGES[j]),r_max=float(EDGES[j+1]),volume=float(total[j,0]),reference_mass=float(total[j,1]),numerical_mass=float(total[j,2]),absolute_density_error=float(total[j,3]),minimum_density=float(minimum[j]),maximum_relative_density_error=float(maximum[j]),worst=worst[j]) for j in range(5)]
    result=dict(t=cp['t'],order=order,input_sha256=hashlib.sha256(data).hexdigest(),bands=rows,
        whole_domain_relative_L1=float(total[:,3].sum()/total[:,1].sum()),
        full_validation=False,scope='Radial diagnostic only; no region removed from full-domain error.')
    atomic('report.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('configuration');p.add_argument('--order',type=int,default=16)
    a=p.parse_args();run(a.source,a.configuration,a.order)
