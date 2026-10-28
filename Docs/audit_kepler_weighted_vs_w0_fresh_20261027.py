"""Compare same-IC W0 and fresh weighted native checkpoint metrics."""

import argparse
import json
from pathlib import Path

from prepare_kepler_sparse_area_restart_20261027 import digest


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('w0',type=Path)
    ap.add_argument('weighted',type=Path)
    ap.add_argument('w0_entropy',type=Path)
    ap.add_argument('weighted_entropy',type=Path)
    ap.add_argument('--allow-w0-failed-check',action='append',default=[],
                    help='Compare while explicitly recording these failed W0 gates')
    ap.add_argument('--allow-weighted-failed-check',action='append',default=[],
                    help='Compare while explicitly recording these failed weighted gates')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    paths=(args.w0,args.weighted,args.w0_entropy,args.weighted_entropy)
    a,b,ea,eb=(json.loads(p.read_text()) for p in paths)
    w0_failed={key for key,value in a['checks'].items() if not value}
    weighted_failed={key for key,value in b['checks'].items() if not value}
    if (w0_failed!=set(args.allow_w0_failed_check) or
            weighted_failed!=set(args.allow_weighted_failed_check) or
            a['initial_sha256']!=b['initial_sha256'] or
            abs(a['time']-b['time'])>1e-6 or
            ea['endpoint_sha256']!=a['endpoint_sha256'] or
            eb['endpoint_sha256']!=b['endpoint_sha256'] or
            ea.get('native_failed_checks',[])!=sorted(w0_failed) or
            eb.get('native_failed_checks',[])!=sorted(weighted_failed) or
            ea['initial_sha256']!=a['initial_sha256'] or
            eb['initial_sha256']!=b['initial_sha256']):
        raise ValueError('Native physics, common initial state, time or entropy lineage mismatch')
    metrics={}
    for field in ('disk_density_abs_mass_weighted','disk_entropy_log_abs_mass_weighted',
                  'disk_position_error_abs_mass_weighted','total_energy_relative'):
        metrics[field]={'w0':a[field],'weighted':b[field],
                        'weighted_over_w0':b[field]/a[field] if a[field] else None}
    metrics['disk_entropy_log_signed_mass_weighted']={
        'w0':ea['whole_disk']['entropy_log_signed_mass_weighted'],
        'weighted':eb['whole_disk']['entropy_log_signed_mass_weighted']}
    metrics['inner_entropy_log_signed_mass_weighted']={
        'w0':ea['radial_bands'][0]['entropy_log_signed_mass_weighted'],
        'weighted':eb['radial_bands'][0]['entropy_log_signed_mass_weighted']}
    result={'initial_sha256':a['initial_sha256'],'w0_time':a['time'],'weighted_time':b['time'],
            'time_difference':b['time']-a['time'],
            'inputs':{str(p.resolve()):digest(p) for p in paths},'metrics':metrics,
            'w0_failed_checks':sorted(w0_failed),
            'weighted_failed_checks':sorted(weighted_failed),
            'checks':{'w0_native_gate_passed':not w0_failed,
                      'w0_failures_match_expected':True,
                      'weighted_native_gate_passed':not weighted_failed,
                      'weighted_failures_match_expected':True,
                      'common_initial_hash':True,'checkpoint_times_within_1e-6':True,
                      'entropy_lineage':True},
            'scope':'Same jittered initial snapshot and almost same checkpoint time. '
                    'Different weighted and W0 binaries and trajectory/restart histories; '
                    'cell density differs with the power geometry. This is not evidence '
                    'that fixed-mass weighted face transport is physically conservative, '
                    'that entropy state differences are irreversible heating, or that '
                    't>6 or GIZMO t=600 stability is established.'}
    with args.output.open('x') as out:
        json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
