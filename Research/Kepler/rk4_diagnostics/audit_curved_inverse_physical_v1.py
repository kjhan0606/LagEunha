"""Exact-decoded fixed-state ledger for inverse target integration.

No tolerance/normalizer changes and no claim about the unassembled targets.
The possible-incidence share is reported separately from the full physical
target tolerance. All nonnegative errors are summed without cancellation.
"""
import csv
import json
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
from audit_curved_target_budget_v1 import r80, verify_regions, power2

ROOT = Path(__file__).resolve().parent


def rbin(s, bits=113):
    x = F(s)
    if not x:
        return x
    sign = -1 if x < 0 else 1
    x = abs(x)
    e = x.numerator.bit_length() - x.denominator.bit_length()
    if x < power2(e):
        e -= 1
    scaled = x * power2(bits - 1 - e)
    n, rem = divmod(scaled.numerator, scaled.denominator)
    if rem * 2 > scaled.denominator or (rem * 2 == scaled.denominator and n % 2):
        n += 1
    return sign * F(n) * power2(e - bits + 1)


def rounded80(x):
    return round_fraction(x, 64)


def round_fraction(x, bits):
    if not x:
        return x
    sign, x = (-1, -x) if x < 0 else (1, x)
    e = x.numerator.bit_length() - x.denominator.bit_length()
    if x < power2(e):
        e -= 1
    scaled = x * power2(bits - 1 - e)
    n, rem = divmod(scaled.numerator, scaled.denominator)
    if rem * 2 > scaled.denominator or (rem * 2 == scaled.denominator and n % 2):
        n += 1
    return sign * F(n) * power2(e - bits + 1)


def read(folder, name):
    with (ROOT / folder / name).open() as handle:
        return list(csv.DictReader(handle))


def dec(x):
    return D(x.numerator) / D(x.denominator)


def string(x):
    return str(dec(x))


def main():
    report = {'scope': 'same225/137 sources; complete local bounds, no all-source target assembly',
              'normalizer_or_physical_tolerance_changed': False,
              'full_cold_transfer_proved': False, 'cases': {}}
    with localcontext() as ctx:
        ctx.prec = 85
        for case in ('225', '137'):
            integral = 'curved_target_integral_v2_' + case
            analytic = 'curved_target_analytic_v1_' + case
            source_glue = 'curved_target_source_glue_v2_' + case
            contracts = {}
            for folder in (integral, analytic, source_glue):
                if (ROOT / folder / 'native_exit_code.txt').read_text().strip() != '0':
                    raise RuntimeError('incomplete run: ' + folder)
                contracts[folder] = json.loads((ROOT / folder / 'contract.json').read_text())
            check, proofs = verify_regions('curved_target_lower_v1_' + case)
            if not check['area_lower_and_image_disjointness_checks']:
                raise RuntimeError('normalizer certificate mismatch')
            counts = {int(r['target']): int(r['possible_element_contributions'])
                      for r in read('curved_error_ledger_v3_' + case, 'incidence_count.csv')}
            flow = read(integral, 'physical_target_moments.csv')
            source = int(flow[0]['source'])
            baseline = {r['moment']: r for r in read(source_glue, 'source_baseline_moments.csv')}
            remainder = {(int(r['target']), r['moment']): r for r in read(analytic, 'analytic_remainders.csv')}
            glue = {(int(r['target']), r['moment']): r for r in read(source_glue, 'target_gluing_bounds.csv')}
            rows, by = [], {}
            for row in flow:
                target, moment = int(row['target']), row['moment']
                analytic_row = remainder[(target, moment)]
                assert analytic_row['all_sectors_complete'] == '1'
                by[(target, moment)] = {
                    'nominal': r80(row['native']),
                    'field_geometry_node_sum': r80(row['field_geometry_node_sum_bound']),
                    'target_quadrature_remainder': rbin(analytic_row['analytic_remainder']),
                    'source_baseline': F(0), 'composition_rounding': F(0),
                    'gluing': F(0), 'represented_primitives': 1,
                }
            moments = list(baseline)
            receivers = sorted(set(int(r['target']) for r in flow))
            for moment in moments:
                b = baseline[moment]
                assert b['analytic_remainder_complete'] == '1'
                own = {'nominal': r80(b['native']), 'field_geometry_node_sum': F(0),
                       'target_quadrature_remainder': F(0),
                       'source_baseline': rbin(b['total_error_bound']),
                       'composition_rounding': F(0), 'gluing': F(0),
                       'represented_primitives': 1 + len(receivers)}
                for receiver in receivers:
                    incoming = by[(receiver, moment)]
                    exact = own['nominal'] - incoming['nominal']
                    rounded = rounded80(exact)
                    own['composition_rounding'] += abs(rounded - exact)
                    own['nominal'] = rounded
                    for key in ('field_geometry_node_sum', 'target_quadrature_remainder'):
                        own[key] += incoming[key]
                by[(source, moment)] = own
            for (target, moment), item in sorted(by.items()):
                grow = glue[(target, moment)]
                assert grow['all_touched_edges_complete'] == '1'
                item['gluing'] = rbin(grow['domain_and_field_ambiguity_bound'])
                if moment not in ('K', 'U', 'potential', 'thermal_rate'):
                    continue
                U = r80(proofs[target]['U_lower'])
                tau = U * (F(1, 10**9) if moment == 'thermal_rate' else F(1, 10**10))
                partial = tau * F(item['represented_primitives'], counts[target])
                components = {k: item[k] for k in ('field_geometry_node_sum', 'target_quadrature_remainder',
                                                   'source_baseline', 'composition_rounding', 'gluing')}
                total = sum(components.values(), F(0))
                rows.append({'target': target, 'moment': moment,
                             'native_partial_contribution': string(item['nominal']),
                             'proved_U_lower': string(U),
                             'physical_whole_target_tolerance': string(tau),
                             'fixed_incidence_partial_allocation': string(partial),
                             'possible_primitives': counts[target],
                             'represented_primitives': item['represented_primitives'],
                             'error_components': {k: string(v) for k, v in components.items()},
                             'complete_local_error_sum': string(total),
                             'local_sum_over_partial_allocation': string(total / partial),
                             'local_sum_over_whole_target_tolerance': string(total / tau),
                             'dominant_component': max(components, key=components.get),
                             'partial_work_allocation_pass': total <= partial,
                             'all_source_target_assembly_pass': False})
            old = {(int(r['target']), r['moment']): r for r in read('curved_fold_transfer_v5_' + case,
                                                                   'target_contribution_bounds.csv')}
            comparison = []
            for receiver in receivers:
                for moment in ('K', 'U', 'potential', 'thermal_rate'):
                    item = by[(receiver, moment)]
                    delta = r80(old[(receiver, moment)]['nominal']) - item['nominal']
                    # This NEW continuous integral enclosure alone validates
                    # the old nominal's deviation; no old sensitivity proof.
                    error = item['field_geometry_node_sum'] + item['target_quadrature_remainder']
                    comparison.append({'target': receiver, 'moment': moment,
                                       'old_minus_new_native': string(delta),
                                       'old_deviation_upper_from_new_integral_enclosure': string(abs(delta) + error),
                                       'new_integral_radius_before_gluing': string(error)})
            worst = max(rows, key=lambda r: D(r['local_sum_over_partial_allocation']))
            worst_whole = max(rows, key=lambda r: D(r['local_sum_over_whole_target_tolerance']))
            report['cases'][case] = {'contracts': contracts, 'normalizer_check': check,
                                     'rows': rows, 'worst_partial': worst,
                                     'worst_whole': worst_whole, 'old_new_comparison': comparison,
                                     'local_terms_complete': True,
                                     'full_target_assembly_complete': False,
                                     'exact_global_C0_proved': False,
                                     'gluing_is_an_explicit_unassigned_strip_error_envelope': True}
            print(case, 'worst partial', worst['target'], worst['moment'],
                  worst['local_sum_over_partial_allocation'], 'dominant', worst['dominant_component'],
                  'worst whole', worst_whole['local_sum_over_whole_target_tolerance'])
    with (ROOT / 'CURVED_INVERSE_PHYSICAL_LEDGER_20261021.json').open('x') as handle:
        json.dump(report, handle, indent=2)


if __name__ == '__main__':
    main()
