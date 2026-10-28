"""Compare native v3/v4 fresh weighted one-step physical fields."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


N=55351
FIELDS=('x','y','vx','vy','mass','volume','ie','den','pressure','w2')
V3='7358c54996b820a7b45bd53bcf6f788f125754b3158845d2b4265afead570f38'
V4='3120aea6dad944b144e04d4597cb3f644abc6c803ff2318b5862b6bf8eec0adc'


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('reference',type=Path)
    ap.add_argument('candidate',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    a=json.loads(args.reference.read_text());b=json.loads(args.candidate.read_text())
    if (not all(a['checks'].values()) or not all(b['checks'].values())
            or a['binary_sha256']!=V3 or b['binary_sha256']!=V4
            or abs(a['time']-.001)>1e-12 or abs(b['time']-.001)>1e-12
            or a['reference_sha256']!=b['reference_sha256']):
        raise ValueError('Native one-step receipt, binary or initial lineage mismatch')
    arrays={}
    for label,d in (('v3',a),('v4',b)):
        for key in ('initial','endpoint'):
            p=Path(d[key])
            if digest(p)!=d[key+'_sha256']:
                raise ValueError(f'Changed {label} {key} snapshot')
            arrays[(label,key)],_=ordered(p,N)
    comparison={}
    for phase in ('initial','endpoint'):
        comparison[phase]={}
        for key in FIELDS:
            x=arrays[('v3',phase)][key];y=arrays[('v4',phase)][key]
            comparison[phase][key]={'bitwise_equal':bool(np.array_equal(x,y)),
                'max_abs_difference':float(np.max(np.abs(x.astype(float)-y.astype(float))))}
    checks={'all_native_gates':True,
            'initial_physical_fields_bitwise':all(v['bitwise_equal'] for v in comparison['initial'].values()),
            'endpoint_physical_fields_bitwise':all(v['bitwise_equal'] for v in comparison['endpoint'].values())}
    result={'reference_receipt':str(args.reference.resolve()),'reference_receipt_sha256':digest(args.reference),
            'candidate_receipt':str(args.candidate.resolve()),'candidate_receipt_sha256':digest(args.candidate),
            'time':.001,'v3_binary_sha256':V3,'v4_binary_sha256':V4,
            'v3_endpoint_sha256':a['endpoint_sha256'],'v4_endpoint_sha256':b['endpoint_sha256'],
            'fields':comparison,'checks':checks,
            'scope':'One native fresh weighted jitter step with nonzero projection gradients; '
                    'exact-zero gradient is covered by the host projection unit test, not this run.'}
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'checks':checks},indent=2))
    if not all(checks.values()):
        raise SystemExit('v3/v4 native physical fields differ')


if __name__=='__main__':
    main()
