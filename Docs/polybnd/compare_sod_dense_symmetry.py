"""Compare transverse moments at the same time on preserved accepted paths."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import restore
from extended_dense_moments import dense_state


def locate_step(root,t,steps):
    lo,hi=1,steps
    while lo<hi:
        mid=(lo+hi)//2;r=json.loads((root/f'step_{mid}_report.json').read_text())
        # The accepted report stores the endpoint time, unlike record['t'].
        if r['t']<t:lo=mid+1
        else:hi=mid
    return lo


def inspect(root,t):
    root=Path(root);cp=json.loads((root/'checkpoint.json').read_text())
    initial=json.loads((root/'initial_state.json').read_text());mass=np.asarray(initial['mass'])
    if t>cp['t'] or t<initial['t']:raise ValueError('requested time outside accepted history')
    lo=locate_step(root,t,cp['steps'])
    p=root/f'step_{lo}_trial.json';data=p.read_bytes();trial=json.loads(data);record=restore(trial['record'])
    assert trial['commit']['accepted'] and trial['commit']['checkpoint']['steps']==lo
    np.testing.assert_array_equal(mass,trial['commit']['checkpoint']['mass'])
    if 'geometry_lift' in record or record.get('dense_moment_arithmetic')!='longdouble':raise ValueError('unsupported path type')
    s=(t-record['t'])/record['dt']
    if not 0<=s<=1:raise ValueError('selected step does not cover requested time')
    y,L,derivative=dense_state(record['y'],record['L'],record['d0'],record['dm'],record['diag'],mass,record['dt'],s)
    a=y[:-4].reshape(-1,20);vy=a[:,10]/mass;i=int(np.argmax(abs(vy)))
    return dict(root=str(root),t=t,step=lo,s=s,trial_sha256=hashlib.sha256(data).hexdigest(),max_mean_vy=float(abs(vy[i])),worst_site=a[i,:2].astype(float).tolist(),transverse_kinetic=float(.5*a[:,18].sum()),min_transverse_second_moment=float(a[:,18].min()))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--time',type=float,required=True);p.add_argument('roots',nargs='+');a=p.parse_args()
    rows=[inspect(root,a.time) for root in a.roots]
    out=dict(t=a.time,rows=rows,full_validation=False,scope='Equal-time dense saved-state diagnostic, no evolution or modified states.')
    Path('report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
