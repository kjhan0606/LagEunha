"""Check field equality and report matched uninstrumented tessellation timing."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest

N=55351
FIELDS=('x','y','vx','vy','mass','volume','ie','den','pressure','w2')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('root',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    cases={}
    good=True
    for case in ('w0','weighted'):
        rows={}; arrays={}
        for mode in ('baseline','enabled'):
            run=args.root/case/mode
            receipt=json.loads((run/'audit.json').read_text())
            if not all(receipt['checks'].values()):
                raise ValueError(f'{case}/{mode} native gate failed')
            endpoint=Path(receipt['endpoint'])
            parent=Path(receipt.get('restart',receipt.get('parent')))
            parent_sha=receipt.get('restart_sha256',receipt.get('parent_sha256'))
            if digest(endpoint)!=receipt['endpoint_sha256'] or digest(parent)!=parent_sha:
                raise ValueError(f'{case}/{mode} snapshot hash changed')
            arrays[mode],_=ordered(endpoint,N)
            log=(run/'log').read_text(errors='replace')
            if '[GPU_SEED_BIN_CERT]' in log:
                raise ValueError('Timing must not enable diagnostic atomics')
            tess=[float(x) for x in re.findall(r'\[GPU tess\].*time=([\d.]+) ms',log)]
            if len(tess)<40:
                raise ValueError(f'{case}/{mode}: too few tessellation samples')
            hot=np.array(tess[4:])
            timer=re.search(r'\[RUN_TIME\] wall=([\d.]+)',(run/'runtime.log').read_text())
            if not timer: raise ValueError('Missing wall time')
            rows[mode]={'receipt':str((run/'audit.json').resolve()),
                        'parent_sha256':parent_sha,'endpoint_sha256':receipt['endpoint_sha256'],
                        'wall_seconds':float(timer.group(1)),'tess_calls':len(tess),
                        'tess_post_warmup_median_ms':float(np.median(hot)),
                        'tess_post_warmup_p90_ms':float(np.quantile(hot,.9))}
        if rows['baseline']['parent_sha256']!=rows['enabled']['parent_sha256']:
            raise ValueError('Different parent')
        fields={}
        for name in FIELDS:
            x,y=arrays['baseline'][name],arrays['enabled'][name]
            equal=bool(np.array_equal(x,y)); good &= equal
            fields[name]={'bitwise_equal':equal,
                          'max_abs_difference':float(np.max(np.abs(x.astype(float)-y.astype(float))))}
        rows['enabled_vs_baseline_tess_median_ratio']=(
            rows['enabled']['tess_post_warmup_median_ms']/rows['baseline']['tess_post_warmup_median_ms'])
        cases[case]={'runs':rows,'fields':fields}
    result={'cases':cases,'checks':{'all_physics_gates':True,'ten_fields_bitwise_equal':bool(good)},
            'scope':'Same GPU, sequential short native replays, diagnostics off. '
                    'Tessellation timings include host and GPU; warmup and ordering '
                    'effects remain. Neither long-time nor MPI equivalence is established.'}
    with args.output.open('x') as out:
        json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({'checks':result['checks'],
                      'tess_median_ratios':{k:v['runs']['enabled_vs_baseline_tess_median_ratio']
                                            for k,v in cases.items()}},indent=2))
    if not good:raise SystemExit('Candidate changed physical fields')


if __name__=='__main__':main()
