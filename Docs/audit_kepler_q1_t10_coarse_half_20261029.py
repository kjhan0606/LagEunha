"""Provisional coarse/half comparison while the quarter-step control runs."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import mean, GAMMA_MINUS_ONE
from audit_kepler_q1_restart_dt_bundle_20261029 import region_metrics


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def distribution(v,m):
    return {'mass_mean':mean(v,m),'cell_p99':float(np.percentile(v,99)),
            'cell_max':float(np.max(v))}


def audit(root):
    cfg=json.loads((root/'stage.json').read_text());seed=Path(cfg['parent']['endpoint'])
    if sha(seed)!=cfg['parent']['sha256'] or cfg['start_time']!=10. or cfg['target_time']!=10.1:
        raise ValueError('Unexpected accepted restart')
    p0,m0=ordered_global(seed);ids=p0['u4if']&ID_MASK
    initial,_=ordered_global(Path(cfg['initial']),ids);r0=np.hypot(initial['x']-3.,initial['y']-3.)
    masks={'core_r0_lt_0p1':r0<.1,'core_r0_0p1_to0p5':(r0>=.1)&(r0<.5),
           'disk_r0_0p5_to2':(r0>=.5)&(r0<2.)}
    states={};cases={};hashes={str(seed):sha(seed),str(Path(cfg['initial'])):sha(cfg['initial'])}
    for name in ['coarse','half']:
        d=root/name;gate=json.loads((d/'checkpoint_gate.json').read_text())
        complete=json.loads((d/'completed.json').read_text());endpoint=Path(gate['endpoint'])
        if (len(gate['checks'])!=28 or not all(gate['checks'].values()) or gate['time']!=10.1
                or complete['sha256']!=sha(endpoint) or sha(d/seed.name)!=cfg['parent']['sha256']
                or sha(d/'eunha2.exe')!=cfg['binary_sha256']):
            raise ValueError('Case failed same-restart, executable or original physics checks: '+name)
        p,meta=ordered_global(endpoint,ids)
        if meta['t']!=10.1 or not np.array_equal(p0['mass'],p['mass']):raise ValueError('Wrong time or mass')
        rows=[(float(re.search(r't_start=(\S+)',line)[1]),float(re.search(r'\bh=(\S+)',line)[1]))
              for line in (d/'log').read_text(errors='replace').splitlines() if '[WALL_ACCEPT]'in line]
        hs=[h for _,h in rows];cap=cfg['requested_caps'][name]
        if (not rows or rows[0][0]!=10. or abs(sum(rows[-1])-10.1)>1e-12
                or any(abs(t+h-rows[k+1][0])>1e-12 for k,(t,h) in enumerate(rows[:-1]))
                or not all(np.isfinite(h) and h>0 and (cap is None or h<=cap*(1+1e-12)) for h in hs)):
            raise ValueError('Invalid native timestep history or requested cap')
        cases[name]={'endpoint':str(endpoint),'sha256':complete['sha256'],'steps':len(hs),
                     'median_h':float(np.median(hs)),
                     'regions':{label:region_metrics(initial,p0,p,mask,10.1) for label,mask in masks.items()}}
        states[name]=p;hashes[str(endpoint)]=complete['sha256']
    if cases['half']['median_h']/cases['coarse']['median_h']>=.75:
        raise ValueError('Actual timestep levels were not separated')
    a,b=states['coarse'],states['half']
    position=np.hypot(a['x']-b['x'],a['y']-b['y'])/np.maximum(r0,.01175)
    density=np.abs(a['den'].astype(float)/b['den'].astype(float)-1)
    entropy=np.abs(np.log(a['ie'].astype(float)/b['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(a['volume'].astype(float)/b['volume'].astype(float)))
    differences={}
    for label,mask in masks.items():
        m=p0['mass'][mask].astype(float)
        differences[label]={'position_difference_over_max_r0_eps':distribution(position[mask],m),
                            'density_relative_difference':distribution(density[mask],m),
                            'entropy_log_difference':distribution(entropy[mask],m)}
    return {'passed':True,'source_hashes':hashes,'cases':cases,'coarse_versus_half':differences,
            'scope':'Two completed native controls from the same inherited t10 state over 0.1 time unit. Provisional before quarter-step completion; not formal convergence order, inherited-error repair, global physical accuracy, or long-time stability.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=audit(args.root)
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))
