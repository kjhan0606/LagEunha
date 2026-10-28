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
common.RANKS.update({name: 1 for name in ['qualifier_old', 'qualifier_new', 'control40', 'vertex40', 'vertex20']})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qualifier', action='store_true')
    args = parser.parse_args()
    if args.qualifier:
        result = common.compare('qualifier_old', 'qualifier_new', exact=True)
        output = ROOT/'qualifier.json'
    else:
        pairs=[]
        for name, halo in [('vertex40', 40), ('vertex20', 20)]:
            state = common.compare('control40', name, exact=True, same_padding=(halo == 40))
            geometry = [compare_dumps(ROOT/'control40'/f'geom_force_rank000_call{call:03d}.bin',
                                      ROOT/name/f'geom_force_rank000_call{call:03d}.bin') for call in range(1,5)]
            text=(ROOT/name/'log').read_text(errors='replace')
            records=re.findall(r'\[VERTEX_HALO\] recertified=(\d+) incomplete=(\d+)',text)
            pad=halo*.0625
            domain=f'actual_rectangle=[{-pad:.17g},{6+pad:.17g}]x[{-pad:.17g},{6+pad:.17g}] W0_only=1'
            checks={'native_trajectory_parity':state['passed'],
                    'all_four_RK4_geometry_and_force_calls_exact':all(x['passed'] for x in geometry),
                    'force_certificate_used_actual_padding_bounds':'[VERTEX_HALO] force=1 '+domain in text,
                    'final_volume_certificate_used_actual_padding_bounds':'[VERTEX_HALO] final_volume=1 '+domain in text,
                    'all_reported_certificates_complete':bool(records) and all(int(x[1]) == 0 for x in records),
                    'smaller_halo_exercised_new_certificate':halo == 40 or any(int(x[0]) > 0 for x in records)}
            pairs.append({'case':name,'checks':checks,'passed':all(checks.values()),
                          'state_comparison':state,'four_force_calls':geometry,
                          'certified_cells_not_covered_by_old_bound_total':sum(int(x[0]) for x in records)})
        qualifier=json.loads((ROOT/'qualifier.json').read_text())
        result={'scope':'Same-checkpoint native geometry/operator and short trajectory test; not t>6 stability or production promotion',
                'qualifier_passed':qualifier['passed'],'pairs':pairs,
                'passed':bool(qualifier['passed'] and all(x['passed'] for x in pairs)),
                'timing_caveat':'Native wall time includes four read-only geometry/force dumps in each full case. No-dump long-run throughput remains to be checked.',
                'production_promoted':False}
        output=ROOT/'comparison.json'
    result['auditor_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']: raise SystemExit('Vertex halo comparison failed; evidence retained')


if __name__ == '__main__': main()
