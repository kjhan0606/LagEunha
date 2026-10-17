"""Equal-time mass-path refinement from independently cached quadrature nodes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integral(v,h):return h/(3*(len(v)-1))*(v[0]+v[-1]+4*v[1:-1:2].sum(0)+2*v[2:-1:2].sum(0))


def prefix(root,audit,t):
    root=Path(root);audit=Path(audit);ic=json.loads((root/'initial_state.json').read_text());mass=np.asarray(ic['mass'])
    identity=json.loads((audit/'identity.json').read_text())
    if identity['initial_sha256']!=digest(root/'initial_state.json') or identity['configuration_sha256']!=digest(root/'configuration.json'):raise ValueError('audit identity mismatch')
    absolute=np.zeros(len(mass));signed=absolute.copy();quad=absolute.copy();time=ic['t'];step=0
    while time<t-1e-15:
        step+=1;p=root/f'step_{step:05d}.npz';nodes=json.loads((audit/f'step_{step:05d}_nodes.json').read_text())
        if nodes['key']!={'identity':identity,'record_sha256':digest(p)}:raise ValueError('cache mismatch')
        if len(nodes['rows'])!=identity['intervals']+1:raise ValueError('incomplete quadrature')
        with np.load(p,allow_pickle=False) as record:
            if abs(float(record['t'])-time)>1e-15:raise ValueError('time gap')
            np.testing.assert_array_equal(record['raw'][:,0],mass);np.testing.assert_array_equal(record['end_raw'][:,0],mass)
            h=float(record['dt']);time+=h
        rates=np.asarray([r['flux'] for r in nodes['rows']]);inc=integral(rates,h)
        absolute+=integral(abs(rates),h);signed+=inc;quad+=abs(inc-integral(rates[::2],h))
    if abs(time-t)>1e-15:raise ValueError('target must be an accepted endpoint')
    return dict(t=time,steps=step,absolute_mass=float(np.max(absolute/mass)),signed_mass=float(np.max(abs(signed)/mass)),nested_difference=float(np.max(quad/mass)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('coarse');p.add_argument('coarse_audit');p.add_argument('fine');p.add_argument('fine_audit');p.add_argument('--times',nargs='+',type=float,required=True);a=p.parse_args()
    left=json.loads((Path(a.coarse)/'initial_state.json').read_text());right=json.loads((Path(a.fine)/'initial_state.json').read_text())
    for k in ('sites','weights','raw','factors','mass'):np.testing.assert_array_equal(left[k],right[k])
    rows=[]
    for t in a.times:
        c=prefix(a.coarse,a.coarse_audit,t);f=prefix(a.fine,a.fine_audit,t)
        rows.append(dict(t=t,coarse=c,fine=f,absolute_mass_ratio=c['absolute_mass']/f['absolute_mass']))
    report=dict(rows=rows,full_validation=False,scope='Equal-time cached mass-path estimates; physical-field accuracy separate.')
    Path('report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
