"""Native W=0 vertex-certificate qualifier, geometry/force and trajectory gate."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import audit_kepler_padding_reuse_pair_20261028 as common
from audit_kepler_geometry_force_dump_20261029 import compare_dumps

ROOT = Path(__file__).resolve().parent
common.ROOT = ROOT
EXPECTED = {'qualifier_old': (0, 40, 0, 0), 'qualifier_new': (0, 40, 0, 0),
            'control40': (0, 40, 0, 0), 'vertex40': (1, 40, 0, 0),
            'vertex20': (1, 20, 0, 0), 'vertex20_reuse': (1, 20, 1, 0),
            'vertex20_reuse_local_sync': (1, 20, 1, 1)}
common.RANKS.update({name: 1 for name in EXPECTED})


def path_checks(name, text):
    enabled, halo, reuse, alias = EXPECTED[name]
    records=re.findall(r'\[VERTEX_HALO\] recertified=(\d+) incomplete=(\d+)',text)
    pad=halo*.0625
    domain=f'actual_rectangle=[{-pad:.17g},{6+pad:.17g}]x[{-pad:.17g},{6+pad:.17g}] W0_only=1'
    checks={'force_certificate_used_actual_padding_bounds':'[VERTEX_HALO] force=1 '+domain in text,
            'final_volume_certificate_used_actual_padding_bounds':'[VERTEX_HALO] final_volume=1 '+domain in text,
            'all_reported_certificates_complete':bool(records) and all(int(x[1]) == 0 for x in records),
            'smaller_halo_exercised_new_certificate':halo == 40 or any(int(x[0]) > 0 for x in records)}
    if reuse:
        checks['padding_reuse_executed']=bool(re.search(r'\[TESS_PADDING_REUSE\] rank=0 .*reuse=1 refresh_validation=1',text))
    if alias:
        checks['device_alias_executed']='[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh' in text
        checks['no_host_refresh_or_alias_fallback']='[HALO_SYNC]' not in text and 'alias invariants failed' not in text
    return checks, records


def launch_checks(root):
    keys=['GFS_GPU_VERTEX_HALO_CERT','GFS_KEPLER_HALO_CELLS','GFS_REUSE_TESS_PADDING','GFS_GPU_LOCAL_PERIODIC_SYNC']
    contexts=[]; checks={}
    for name, expected in EXPECTED.items():
        path=root/name/'launch_context.json'
        checks[name+'_launch_context_present']=path.is_file()
        if not path.is_file(): continue
        context=json.loads(path.read_text());contexts.append(context)
        checks[name+'_requested_flags_exact']=all(context.get(key)==str(value) for key,value in zip(keys,expected))
    checks['all_seven_cases_recorded']=len(contexts)==len(EXPECTED)
    jobs={context.get('SLURM_JOB_ID','') for context in contexts}
    checks['one_valid_Slurm_allocation']=len(jobs)==1 and all(job.isdecimal() for job in jobs)
    nodes={context.get('SLURM_NODELIST','') for context in contexts}
    checks['one_recorded_node_allocation']=len(nodes)==1 and all(nodes)
    checks['same_visible_GPU_mask']=len({context.get('CUDA_VISIBLE_DEVICES','') for context in contexts})==1
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qualifier', action='store_true')
    args = parser.parse_args()
    if args.qualifier:
        result = common.compare('qualifier_old', 'qualifier_new', exact=True)
        output = ROOT/'qualifier.json'
    else:
        pairs=[]
        for control, name in [('control40', 'vertex40'), ('control40', 'vertex20'),
                              ('vertex20', 'vertex20_reuse'), ('vertex20_reuse', 'vertex20_reuse_local_sync')]:
            state = common.compare(control, name, exact=True, same_padding=(EXPECTED[control][1] == EXPECTED[name][1]))
            geometry = [compare_dumps(ROOT/control/f'geom_force_rank000_call{call:03d}.bin',
                                      ROOT/name/f'geom_force_rank000_call{call:03d}.bin') for call in range(1,5)]
            text=(ROOT/name/'log').read_text(errors='replace')
            checks,records=path_checks(name,text)
            checks.update(native_trajectory_parity=state['passed'],all_four_RK4_geometry_and_force_calls_exact=all(x['passed'] for x in geometry))
            timings={}
            for case in [control,name]:
                log=(ROOT/case/'log').read_text(errors='replace')
                totals=[float(value) for value in re.findall(r'^\[CPU step \d+\].*?total=([0-9.]+)s',log,re.M)]
                timings[case]={'native_CPU_step_mean_after_first_seconds':sum(totals[1:])/len(totals[1:]) if len(totals)>1 else None,
                               'steady_records':max(0,len(totals)-1)}
            pairs.append({'control':control,'case':name,'checks':checks,'passed':all(checks.values()),
                          'steady_internal_timing':timings,
                          'state_comparison':state,'four_force_calls':geometry,
                          'certified_cells_not_covered_by_old_bound_total':sum(int(x[0]) for x in records)})
        qualifier=json.loads((ROOT/'qualifier.json').read_text())
        contexts=launch_checks(ROOT)
        result={'scope':'Same-checkpoint native geometry/operator and short trajectory test; not t>6 stability or production promotion',
                'qualifier_passed':qualifier['passed'],'pairs':pairs,
                'launch_context_checks':contexts,
                'passed':bool(qualifier['passed'] and all(x['passed'] for x in pairs) and all(contexts.values())),
                'timing_caveat':'Wall includes four read-only dumps per full case. Native internal step timer excludes its first instrumented step but does not cover all driver wall cost. No-dump long-run throughput remains to be checked.',
                'production_promoted':False}
        output=ROOT/'comparison.json'
    result['auditor_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']: raise SystemExit('Vertex halo comparison failed; evidence retained')


if __name__ == '__main__': main()
