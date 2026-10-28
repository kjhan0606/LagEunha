"""Resolve signed entropy changes in an audited fixed-mass weighted state."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_entropy_sign_20261027 import EDGES, GM1, N, describe
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('receipt',type=Path)
    ap.add_argument('--allow-failed-check',action='append',default=[],
                    help='Diagnose, but do not hide, these named failed gates')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    receipt=json.loads(args.receipt.read_text())
    required={'finite_positive','fixed_mass_exact','initial_W0_and_evolved_weights',
              'density_and_EOS','zero_injection','total_energy','disk_entropy',
              'independent_power_area','weight_face_guard','no_geometry_abort'}
    failed={key for key,value in receipt['checks'].items() if not value}
    if (len(receipt['checks'])!=21 or not required.issubset(receipt['checks']) or
            failed!=set(args.allow_failed_check)):
        raise ValueError('Weighted parent receipt did not pass the 21 checks')
    initial,endpoint=Path(receipt['initial']),Path(receipt['endpoint'])
    for key,path in (('initial',initial),('endpoint',endpoint)):
        if digest(path)!=receipt[key+'_sha256']:raise ValueError(f'Changed {key}')
    p0,m0=ordered(initial,N);p,m=ordered(endpoint,N)
    if abs(m0['t'])>1e-12 or abs(m['t']-receipt['time'])>1e-10:
        raise ValueError('Checkpoint time changed')
    r=np.hypot(p0['x'].astype(float)-3.,p0['y'].astype(float)-3.)
    selected=(r>=EDGES[0])&(r<EDGES[-1]);radius=r[selected]
    entropy=(np.log(p['ie'][selected].astype(float)/p0['ie'][selected].astype(float))
             +GM1*np.log(p['volume'][selected].astype(float)/p0['volume'][selected].astype(float)))
    mass=p0['mass'][selected].astype(float)
    if not np.all(np.isfinite(entropy)):raise ValueError('Nonfinite entropy')
    whole=describe(entropy,mass)
    if abs(whole['entropy_log_abs_mass_weighted']-receipt['disk_entropy_log_abs_mass_weighted'])>1e-12:
        raise ValueError('Entropy differs from weighted receipt')
    bands=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        mask=(radius>=lo)&(radius<hi)
        bands.append({'r_min':lo,'r_max':hi,**describe(entropy[mask],mass[mask])})
    result={'accepted_receipt':str(args.receipt.resolve()),'initial_sha256':receipt['initial_sha256'],
            'endpoint_sha256':receipt['endpoint_sha256'],'time':receipt['time'],
            'native_failed_checks':sorted(failed),'whole_disk':whole,'radial_bands':bands,
            'scope':'Signed entropy state difference in an experimental fixed-mass weighted '
                    'trajectory; not an irreversible-heating or moving-face flux ledger.'}
    with args.output.open('x') as out:
        json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({'time':result['time'],'whole_disk':whole,'inner_disk':bands[0]},indent=2))


if __name__=='__main__':main()
