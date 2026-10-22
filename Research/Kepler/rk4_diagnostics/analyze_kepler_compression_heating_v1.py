#!/usr/bin/env python3
"""Saved-data diagnostic of the four-step ORIGINAL-rule trajectory only.

No point evaluation, field/force/mass solve or trajectory is executed here.
The Jensen excess is signed: finite-sum/CSV noise is never clamped into heat.
"""
import csv
import decimal
import hashlib
import json
import math
import pathlib
import struct
import time

D = decimal.Decimal
decimal.getcontext().prec = 55
ROOT = pathlib.Path(__file__).resolve().parent
INITIAL = ROOT / 'material_full_initial_ledger_v1_mpi4'
TRAJECTORY = ROOT / 'material_full_cached_continuation_v1_ada_mpi4'
STEP = TRAJECTORY / 'run/continuous_step4'
DEST = ROOT / 'material_full_compression_heating_audit_v1'
EDGES = [0, .15, .5, 2, 8, 100]


def read_json(path):
    return json.loads(path.read_text(), parse_float=D)


def blank():
    return dict(cells=0, U0=D(0), U=D(0), uniform_U=D(0),
                min_signed_excess=None, max_signed_excess=None,
                max_abs_bulk=D(0), max_abs_total=D(0))


def update(group, u0, u, uniform, excess, bulk, total):
    group['cells'] += 1
    group['U0'] += u0
    group['U'] += u
    group['uniform_U'] += uniform
    for key, operator in [('min_signed_excess', min), ('max_signed_excess', max)]:
        group[key] = excess if group[key] is None else operator(group[key], excess)
    group['max_abs_bulk'] = max(group['max_abs_bulk'], abs(bulk))
    group['max_abs_total'] = max(group['max_abs_total'], abs(total))


def summarize(group):
    group = dict(group)
    u0 = group['U0']
    assert u0 > 0
    group['total_U_change_over_U0'] = (group['U'] - u0) / u0
    group['uniform_compression_change_over_U0'] = (group['uniform_U'] - u0) / u0
    group['signed_Jensen_excess_over_U0'] = (group['U'] - group['uniform_U']) / u0
    return group


def main():
    begin = time.monotonic()
    assert (INITIAL / 'exit_code.txt').read_text().strip() == '0'
    assert (TRAJECTORY / 'native_exit_code.txt').read_text().strip() == '0'
    assert read_json(ROOT / 'material_full_cached_audit_v1/contract.json')['audit_pass']
    step = read_json(STEP / 'step.json')
    assert step['finite_step_gates_pass'] and step['time_after'] < D('.001')
    for stage in range(1, 5):
        rate = read_json(STEP / f'stage{stage}/rhs.json')
        assert rate['fixed_positive_points_per_pass'] == 2040256512
        assert rate['RHS_gates_pass']
    DEST.mkdir()
    bins = [blank() for _ in range(len(EDGES) - 1)]
    total = blank()
    seen = set()
    files = [INITIAL / 'exit_code.txt', TRAJECTORY / 'native_exit_code.txt',
             ROOT / 'material_full_cached_audit_v1/contract.json', STEP / 'step.json']
    extreme = {'min': None, 'max': None}
    law_reference = None
    with (DEST / 'cells.csv').open('x', newline='') as output:
        writer = csv.writer(output)
        writer.writerow(['cell', 'radius_mass_centroid', 'U_initial', 'volume_initial',
                         'volume_current', 'U_current', 'U_uniform_compression',
                         'bulk_compression_fraction', 'total_U_fraction',
                         'signed_Jensen_excess_fraction'])
        for rank in range(4):
            thermal = ROOT / 'material_full_initial_v2_h200_mpi4/run' / f'full_generating_thermal_r{rank}.bin'
            with thermal.open('rb') as f:
                header = struct.unpack('<6Q', f.read(48))
                law = struct.unpack('<4d', f.read(32))
            assert header[0] == 0x4d46544845524d31 and header[1] == 1
            assert header[3:] == (16, rank, 1)
            if law_reference is None:
                law_reference = law
            assert law == law_reference and law[0] > 1 and law[1] > 0 and law[3] == 0
            gamma, pressure = D.from_float(law[0]), D.from_float(law[1])
            initial = INITIAL / 'run' / f'full_initial_local_ledger_r{rank}.csv'
            endpoint = STEP / 'endpoint' / f'physical_cells_r{rank}.csv'
            files.extend([thermal, initial, endpoint])
            with initial.open() as f:
                before = {int(row['cell']): D(row['U']) for row in csv.DictReader(f)}
            assert len(before) == header[2]
            rank_seen = set()
            with endpoint.open() as f:
                for row in csv.DictReader(f):
                    cell = int(row['cell'])
                    assert cell in before and cell not in seen
                    seen.add(cell)
                    rank_seen.add(cell)
                    u0, u, volume = before[cell], D(row['U']), D(row['volume'])
                    assert u0 > 0 and u > 0 and volume > 0
                    volume0 = u0 * (gamma - 1) / pressure
                    ratio = volume / volume0
                    uniform_factor = ((1 - gamma) * ratio.ln()).exp()
                    uniform = u0 * uniform_factor
                    bulk, change = uniform_factor - 1, u / u0 - 1
                    excess = change - bulk
                    radius = math.hypot(float(row['x_mass']), float(row['y_mass']))
                    bucket = next(i for i in range(len(bins)) if EDGES[i] <= radius < EDGES[i + 1])
                    update(total, u0, u, uniform, excess, bulk, change)
                    update(bins[bucket], u0, u, uniform, excess, bulk, change)
                    writer.writerow([cell, repr(radius), u0, volume0, volume, u, uniform, bulk, change, excess])
                    record = dict(cell=cell, radius=radius, signed_excess=excess, bulk=bulk, total=change)
                    if extreme['min'] is None or excess < extreme['min']['signed_excess']:
                        extreme['min'] = record
                    if extreme['max'] is None or excess > extreme['max']['signed_excess']:
                        extreme['max'] = record
            assert rank_seen == set(before)
    assert total['cells'] == 71080
    report = dict(scope='saved original-rule step4 only; no physical evolution or force/point replay',
                  accepted_time=step['time_after'], cells=total['cells'],
                  law_gamma_exact_double=D.from_float(law_reference[0]),
                  law_pressure_exact_double=D.from_float(law_reference[1]),
                  global_values=summarize(total), radial_edges=EDGES,
                  radial_values=[summarize(g) for g in bins], extreme_cells=extreme,
                  Jensen_identity='U_min=U0*(V/V0)^(1-gamma); signed_excess=(U-U_min)/U0',
                  initial_context='fresh constant P0, same generating entropy and original positive rule',
                  signs_near_finite_sum_or_CSV_roundoff_unresolved=True,
                  negative_excess_clamped=False, pointwise_variance_measured=False,
                  reduced_rule_future_states_need_their_own_initial_measure_context=True,
                  artificial_heating_cause_proven=False, t6_t8_t10_validated=False,
                  seconds=time.monotonic() - begin)
    (DEST / 'contract.json').write_text(json.dumps(report, default=str, indent=2) + '\n')
    with (DEST / 'input.sha256').open('x') as out:
        for path in files:
            out.write(hashlib.sha256(path.read_bytes()).hexdigest() + '  ' + str(path) + '\n')
    (DEST / 'source.sha256').write_text(hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest() + '  ' + str(pathlib.Path(__file__).resolve()) + '\n')
    print(json.dumps(report, default=str, indent=2))


if __name__ == '__main__':
    main()
