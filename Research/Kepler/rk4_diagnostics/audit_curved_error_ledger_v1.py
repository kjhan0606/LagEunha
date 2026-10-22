"""Summarize the fixed error-attribution bundle, never accept a transfer.

Only completed immutable native runs are inputs. Decimal totals preserve
tiny cold budgets; an estimator is kept separate from proven error terms.
"""
import csv
import json
from decimal import Decimal as D, localcontext
from pathlib import Path
from audit_curved_target_budget_v1 import verify_regions

ROOT = Path(__file__).resolve().parent


def read(folder, name):
    with (ROOT / folder / name).open() as handle:
        return list(csv.DictReader(handle))


def main():
    report = {
        'scope': 'same225/137/lens saved positive packets; attribution only',
        'source_of_nominal_state_changed': False,
        'full_cold_transfer_proved': False,
        'missing_certificates': [
            'source-reference baseline arithmetic and rigorous integration remainder',
            'quadrature integration remainder (retained sensitivity is an estimator)',
            'whole shared-edge gap to target-local domain and field influence',
            'all possible contributors of one assembled physical target',
        ],
        'dependency_paths': [
            {'path': 'regular labels and roots', 'reference': 'compiled midpoint representative',
             'charged': 'original conversion + actual scalar graph + root guard',
             'native_position_as_boundary_shift': False},
            {'path': 'rounded power-plane scalars',
             'charged': 'existing exact-plane to rounded-plane conversion',
             'charged_again_in_actual_ray_graph': False},
            {'path': 'event unknown owner', 'reference': 'native position',
             'original_event_protection_retained': True},
            {'path': 'integrand', 'reference': 'original stored-coefficient mathematical basis',
             'charged': 'native fields + rho + all12 moments including K,U,Udot,potential'},
            {'path': 'pre-save rho and weight',
             'charged': 'binary80 interval within half a binary64 ulp; all original operations'},
            {'path': 'original accumulation',
             'charged': 'binary80 gamma_(target points + accepted panels +4) times absolute term sum'},
        ],
        'cases': {},
    }
    components = {
        'field_forward_bound': 'field',
        'serialization_bridge_bound': 'serialization',
        'old_binary80_sum_bound': 'sum',
        'updated_regular_boundary_bound': 'regular_boundary',
        'event_bound_retained': 'event',
        'reference_reconstruction_radius': 'reference_oracle',
    }
    with localcontext() as ctx:
        ctx.prec = 80
        for case in ('225', '137', 'lens'):
            folder = 'curved_error_ledger_v3_' + case
            if (ROOT / folder / 'native_exit_code.txt').read_text().strip() != '0':
                raise RuntimeError('incomplete native attribution: ' + case)
            point = json.loads((ROOT / ('curved_original_oracle_v2_' + case) / 'contract.json').read_text())
            if not point['point_forward_model_pass']:
                raise RuntimeError('independent point model failed: ' + case)
            normalizers, proofs = {}, []
            for proof in ['curved_target_lower_v1_' + case] + (['curved_target_lower_v2_lens'] if case == 'lens' else []):
                check, bounds = verify_regions(proof)
                if not check['area_lower_and_image_disjointness_checks']:
                    raise RuntimeError('positive normalizer proof mismatch')
                proofs.append(check)
                for target, row in bounds.items():
                    if D(row['U_lower']) > D(normalizers.get(target, {}).get('U_lower', '-1')):
                        normalizers[target] = {**row, 'proof_folder': proof}
            counts = {int(r['target']): int(r['possible_element_contributions'])
                      for r in read(folder, 'incidence_count.csv')}
            data = read(folder, 'moment_ledger.csv')
            flow = {}
            for row in data:
                if row['moment'] not in ('K', 'U', 'potential', 'thermal_rate'):
                    continue
                source, receiver = int(row['source']), int(row['target'])
                for affected in (source, receiver):
                    key = (affected, row['moment'])
                    target = flow.setdefault(key, {'components': {k: D(0) for k in components.values()},
                                                  'quad': D(0), 'observed': D(0), 'primitive_contributions': 0})
                    for col, label in components.items():
                        target['components'][label] += D(row[col])
                    target['quad'] += D(row['quadrature_sensitivity_estimate'])
                    target['observed'] += abs(D(row['observed_old_nominal_minus_reference_q128']))
                    target['primitive_contributions'] += 1
            rows = []
            for (target, moment), raw in sorted(flow.items()):
                lower = D(normalizers[target]['U_lower'])
                assert lower > 0 and 0 < raw['primitive_contributions'] <= counts[target]
                tau = lower * (D('1e-9') if moment == 'thermal_rate' else D('1e-10'))
                partial = tau * D(raw['primitive_contributions']) / D(counts[target])
                known = sum(raw['components'].values(), D(0))
                dominant = max(raw['components'], key=raw['components'].get)
                rows.append({
                    'target': target, 'moment': moment, 'U_lower': str(lower),
                    'proof_folder': normalizers[target]['proof_folder'],
                    'possible_primitive_contributions': counts[target],
                    'represented_primitive_contributions': raw['primitive_contributions'],
                    'physical_target_tolerance': str(tau),
                    'fixed_incidence_partial_allocation': str(partial),
                    'known_error_model_components': {k: str(v) for k, v in raw['components'].items()},
                    'known_error_model_sum': str(known),
                    'known_sum_over_partial_allocation': str(known / partial),
                    'known_sum_over_whole_target_tolerance': str(known / tau),
                    'dominant_known_term': dominant,
                    'dominant_term_over_partial_allocation': str(raw['components'][dominant] / partial),
                    'quadrature_sensitivity_estimate': str(raw['quad']),
                    'quadrature_estimate_over_partial_allocation': str(raw['quad'] / partial),
                    'observed_same_packet_native_minus_original_reference_abs_sum': str(raw['observed']),
                    'observed_same_packet_difference_over_partial_allocation': str(raw['observed'] / partial),
                    'observed_same_packet_difference_is_physical_transfer_error': False,
                    'source_and_gluing_bounds': None,
                    'acceptance': False,
                })
            worst = max(rows, key=lambda row: D(row['known_sum_over_partial_allocation']))
            report['cases'][case] = {
                'point_oracle': point,
                'native_ledger': json.loads((ROOT / folder / 'contract.json').read_text()),
                'normalizer_proofs': proofs,
                'rows': rows,
                'worst_partial_row': worst,
                'maximum_known_over_whole_target_tolerance': str(max(D(row['known_sum_over_whole_target_tolerance']) for row in rows)),
                'maximum_observed_same_packet_over_allocation': str(max(D(row['observed_same_packet_difference_over_partial_allocation']) for row in rows)),
                'full_cold_transfer_proved': False,
            }
            print(case, 'worst target', worst['target'], worst['moment'],
                  'known/allocation', worst['known_sum_over_partial_allocation'],
                  'dominant', worst['dominant_known_term'],
                  'max-known/full-tolerance', report['cases'][case]['maximum_known_over_whole_target_tolerance'],
                  'observed/allocation-max', report['cases'][case]['maximum_observed_same_packet_over_allocation'])
    with (ROOT / 'CURVED_ORIGINAL_ERROR_LEDGER_20261021.json').open('x') as handle:
        json.dump(report, handle, indent=2)


if __name__ == '__main__':
    main()
