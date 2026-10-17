"""Fixed-geometry source-flow convergence of the current enriched gravity kick."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from audit_kepler_energy_step import configure,phi,atomic
from potential_moment_coupled_2d import provider
from gravity_source_2d import apply
from kepler_disk_reference import polar_rule,frozen,CENTER
from kepler_spline_gravity import acceleration


def run(seed,initial_potential,duration=.001,scheme='finite'):
    if scheme=='midpoint':
        from midpoint_gravity_source import apply as source_apply
    elif scheme=='finite':source_apply=apply
    else:raise ValueError('unknown source scheme')
    root=Path(seed);data={n:(root/n).read_bytes() for n in ('source.json','configuration.json')}
    data['initial_potential.json']=Path(initial_potential).read_bytes()
    identity={k:hashlib.sha256(v).hexdigest() for k,v in data.items()};identity['duration']=duration
    identity['scheme']=scheme
    identity['source_manifest_sha256']=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest()
    if Path('identity.json').exists():
        if json.loads(Path('identity.json').read_text())!=identity:raise ValueError('source audit identity changed')
    else:atomic('identity.json',identity)
    cp=json.loads(data['source.json']);configure(json.loads(data['configuration.json']))
    y=np.asarray(cp['state']);mass=np.asarray(cp['mass']);L=[np.asarray(a) for a in cp['factors']];U=np.asarray(json.loads(data['initial_potential.json']))
    mask=frozen(cp['t'],y[:-4].reshape(-1,20)[:,:2]);get=provider(U,phi,lambda p:polar_rule(p,16,radii=(.25,.5)),mask)
    outputs=[]
    for n in (1,2,4,8,16):
        file=Path(f'subdivision_{n}.json')
        if file.exists():result=json.loads(file.read_text())
        else:
            state=y.copy();factors=[a.copy() for a in L]
            for j in range(n):state,factors,diag=source_apply(cp['t'],state,mass,factors,lambda x:acceleration(x,center=CENTER,rc=.5),duration/n,mask,get)
            result=dict(n=n,state=state,factors=factors);atomic(file,result)
        outputs.append(np.asarray(result['state']));print('gravity subdivision',n,flush=True)
    pairs=[]
    for n,a,b in zip((1,2,4,8),outputs[:-1],outputs[1:]):
        diff=(a-b)[:-4].reshape(-1,20)
        pairs.append(dict(coarse_n=n,max_state_difference=float(np.max(abs(diff))),
            sum_active_energy_difference=float(diff[~mask,19].sum()),max_energy_difference=float(np.max(abs(diff[~mask,19])))))
    report=dict(identity=identity,pairs=pairs,full_validation=False,scope='gravity source flow only; fixed density/geometry and original moments')
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');p.add_argument('initial_potential');p.add_argument('--scheme',choices=('finite','midpoint'),default='finite');a=p.parse_args();run(a.seed,a.initial_potential,scheme=a.scheme)
