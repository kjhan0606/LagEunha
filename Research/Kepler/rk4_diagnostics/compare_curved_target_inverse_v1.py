"""Compare preserved polar moments with the certified new physical target.

This bounds the old moment's actual deviation using the NEW integral
enclosure. It does not turn the old rule difference into a remainder proof.
"""
import csv
import json
from decimal import Decimal as D, localcontext
from pathlib import Path
from audit_curved_target_budget_v1 import r80, dec

ROOT = Path(__file__).resolve().parent


def read(path):
    with (ROOT / path).open() as handle:
        return list(csv.DictReader(handle))


def main():
    folder = ROOT / 'curved_target_inverse_v2_lens'
    if (folder / 'native_exit_code.txt').read_text().strip() != '0':
        raise RuntimeError('new inverse target is not a completed pass')
    old = {r['moment']: r for r in read('curved_fold_transfer_v5_lens/target_contribution_bounds.csv') if r['target'] == '101'}
    report = {'scope': 'same original tiny target101, new validated continuous-integral enclosure',
              'old_packet_preserved': True, 'target_U_lower_changed': False,
              'full_source_or_all_target_transfer': False, 'rows': []}
    with localcontext() as ctx:
        ctx.prec = 80
        for row in read('curved_target_inverse_v2_lens/target_moments.csv'):
            name = row['moment']
            difference = r80(old[name]['nominal']) - r80(row['native'])
            radius = r80(row['total_error_bound'])
            observed = abs(difference)
            lower, upper = max(0, observed - radius), observed + radius
            tolerance = r80(row['unchanged_tolerance'])
            report['rows'].append({
                'moment': name,
                'signed_old_minus_new_native': str(dec(difference)),
                'new_continuous_integral_radius': str(dec(radius)),
                'old_deviation_from_continuous_integral_lower': str(dec(lower)),
                'old_deviation_from_continuous_integral_upper': str(dec(upper)),
                'old_deviation_upper_divided_by_tolerance': str(dec(upper / tolerance)),
                'new_error_bound_divided_by_tolerance': row['budget_ratio'],
                'old_rule_sensitivity_estimate': old[name]['quad_estimate'],
                'old_sensitivity_is_not_a_remainder_proof': True,
            })
        worst = max(report['rows'], key=lambda r: D(r['old_deviation_upper_divided_by_tolerance']))
        report['max_old_continuous_error_over_tolerance'] = worst
        print('old max continuous deviation enclosure / tolerance:', worst['moment'],
              worst['old_deviation_upper_divided_by_tolerance'])
    with (ROOT / 'CURVED_TARGET_INVERSE_COMPARISON_20261021.json').open('x') as handle:
        json.dump(report, handle, indent=2)


if __name__ == '__main__':
    main()
