"""Exact-decoded per-source baseline/outflow/complement conservation.
This is not conservation of the unassembled24-source-cell state.
"""
import json
from decimal import localcontext
from fractions import Fraction as F
from audit_curved_inverse_physical_v1 import ROOT, read, r80, rounded80, string


def main():
    result = {'scope': 'same225/137 per-source represented12 moments',
              'whole24_cell_state_conservation_verified': False, 'cases': {}}
    with localcontext() as ctx:
        ctx.prec = 85
        for case in ('225', '137'):
            baseline = {r['moment']: r for r in read('curved_target_source_glue_v2_' + case,
                                                     'source_baseline_moments.csv')}
            flows = read('curved_target_integral_v2_' + case, 'physical_target_moments.csv')
            rows = []
            for moment, base in baseline.items():
                reference = r80(base['native'])
                own, rounding = reference, F(0)
                incoming = []
                for row in sorted((r for r in flows if r['moment'] == moment), key=lambda r: int(r['target'])):
                    value = r80(row['native'])
                    exact = own - value
                    n = rounded80(exact)
                    rounding += abs(n - exact)
                    own = n
                    incoming.append(value)
                residual = own + sum(incoming, F(0)) - reference
                assert abs(residual) <= rounding
                rows.append({'moment': moment, 'source_baseline': string(reference),
                             'positive_outgoing_sum_signed_moment': string(sum(incoming, F(0))),
                             'source_own_complement': string(own),
                             'source_partition_residual': string(residual),
                             'actual_binary80_composition_rounding_bound': string(rounding),
                             'conservation_within_actual_rounding': True})
            assert next(r for r in rows if r['moment'] == 'mass')['source_own_complement'][0] != '-'
            result['cases'][case] = rows
            print(case, 'all12 source moments conserve within actual binary80 composition rounding')
    with (ROOT / 'CURVED_INVERSE_SOURCE_CONSERVATION_20261021.json').open('x') as handle:
        json.dump(result, handle, indent=2)


if __name__ == '__main__':
    main()
