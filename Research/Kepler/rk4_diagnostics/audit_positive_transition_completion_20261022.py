#!/usr/bin/env python3
"""Read completed native outputs; never evolve, repair, or rewrite a state."""
from pathlib import Path
import csv
import datetime
import hashlib
import json
import math
import struct

ROOT = Path(__file__).resolve().parent
RUN = ROOT / 'material_full_positive_transition_v2_ada_mpi4'
STEP = RUN / 'run/continuous_step5'
DEST = ROOT / 'material_full_positive_transition_completion_audit_v1'
assert ROOT == Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics')
DEST.mkdir(exist_ok=False)
evidence = {}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(chunk)
    evidence[str(path)] = h.hexdigest()
    return h.hexdigest()


def document(path):
    digest(path)
    return json.loads(path.read_text())


native = RUN / 'native_exit_code.txt'
assert native.read_text().strip() == '0'
digest(native)
dispatch = ROOT / 'material_full_positive_dispatch_v1'
for name in ['exit_code.txt', 'ada_runner_exit_code.txt']:
    f = dispatch / name
    assert f.read_text().strip() == '0'
    digest(f)
terminal = document(dispatch / 'state.json')
assert terminal['state'] == 'actual_paired_step_complete'
integrity = {}
for name, count in [
    ('final_source_integrity.txt', 68), ('final_binary_integrity.txt', 2),
    ('final_accepted_input_integrity.txt', 7), ('final_prepared_k1_integrity.txt', 6),
    ('final_positive_rule_integrity.txt', 9), ('immutable_payload_integrity.txt', 29),
]:
    f = RUN / name
    lines = f.read_text().splitlines()
    assert len(lines) == count and all(line.endswith(': OK') for line in lines)
    integrity[name] = count
    digest(f)
contract = document(RUN / 'contract.json')
result = document(RUN / 'run/full_continuous_result.json')
step = document(STEP / 'step.json')
cost = document(STEP / 'selected_rule_cost.json')
restart = document(STEP / 'numerical_restart_roundtrip.json')
physical = document(STEP / 'endpoint/physical_state.json')
assert result['new_accepted_steps'] == 1 and result['total_accepted_steps'] == 5
assert result['same_accepted_trajectory'] and result['no_replayed_accepted_interval']
assert step['time_before'] == result['time_before'] and step['time_after'] == result['accepted_time']
assert math.isclose(step['time_after'] - step['time_before'], step['h'], rel_tol=1e-14)
assert step['finite_step_gates_pass']
for name, bound in [('energy_defect_over_U', 1e-8), ('thermal_increment_defect_over_U', 1e-8),
                    ('angular_defect_over_2K', 1e-11), ('mass_relative_difference', 1e-14),
                    ('mass_vs_point_K_error_over_U', 1e-8)]:
    assert abs(step[name]) < bound
assert restart['version'] == 3 and not restart['accepted_interval_replayed']
for name in ['actual_accepted_checkpoint_reloaded', 'coefficients_hash_match',
             'physical_and_numerical_ledgers_match', 'exact_rule_generation_and_active_mode_match']:
    assert restart[name]
assert physical['time'] == result['accepted_time']
assert all(physical[name] > 0 for name in ['min_J', 'min_jhat', 'min_U_point', 'physical_dt_bound'])
assert physical['no_state_or_heat_correction'] and not physical['global_map_injectivity_proven']

stages = []
cell_ids = None
for number in range(1, 5):
    base = STEP / f'stage{number}'
    comparison = document(base / 'quadrature_difference.json')
    assert comparison['same_stage_quadrature_gate_pass'] and comparison['same_positions_velocities']
    assert comparison['no_physical_state_correction']
    assert comparison['high_k1_force_and_solve_reused'] == (number == 1)
    for quantity in ['pressure', 'density', 'velocity']:
        assert 0 <= comparison[f'cumulative_{quantity}_error'] <= comparison['error_allowance_at_step_end']
    rhs = document(base / 'low/rhs.json')
    assert rhs['RHS_gates_pass'] and rhs['GPU_mass_acceleration_implemented']
    assert max(rhs['mass_true_residual']) < 1e-16
    ids = set()
    max_cell_energy = max_cell_angular = 0.0
    for rank in range(4):
        f = base / 'low' / f'cells_r{rank}.csv'
        with f.open() as stream:
            for row in csv.DictReader(stream):
                cell = int(row['cell'])
                assert cell not in ids
                ids.add(cell)
                assert all(math.isfinite(float(v)) for key, v in row.items() if key != 'cell')
                assert all(float(row[key]) > 0 for key in ['volume', 'density', 'mean_pressure', 'U', 'K', 'inertia'])
                max_cell_energy = max(max_cell_energy, float(row['energy_rate_error_over_U']))
                max_cell_angular = max(max_cell_angular, float(row['angular_rate_error_over_2K']))
        digest(f)
        for mode in ['low', 'high']:
            f = base / mode / f'passes_r{rank}.tsv'
            with f.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            for pass_id in ({0, 1} if mode == 'high' and number == 1 else {0, 1, 2}):
                assert max(int(row['triangles_done']) for row in rows if int(row['pass']) == pass_id) == 284608
            digest(f)
    assert len(ids) == 71080
    assert max_cell_energy < 1e-8 and max_cell_angular < 1e-11
    assert math.isclose(max_cell_energy, rhs['every_cell_energy_rate_error_over_U'], rel_tol=1e-14)
    assert math.isclose(max_cell_angular, rhs['every_cell_angular_rate_error_over_2K'], rel_tol=1e-14)
    if cell_ids is None:
        cell_ids = ids
    assert cell_ids == ids
    stages.append({'stage': number, 'cells_checked': len(ids), 'quadrature': comparison, 'rhs': rhs})

endpoint = document(STEP / 'endpoint/quadrature_difference.json')
assert endpoint['endpoint_quadrature_gate_pass'] and endpoint['same_positions_velocities']
assert not endpoint['high_force_or_solve_repeated'] and not endpoint['transition_difference_is_physical_heat']
assert endpoint['no_physical_state_correction'] and not endpoint['global_map_injectivity_proven']
ids = set()
for rank in range(4):
    f = STEP / 'endpoint' / f'physical_cells_r{rank}.csv'
    with f.open() as stream:
        for row in csv.DictReader(stream):
            cell = int(row['cell'])
            assert cell not in ids
            ids.add(cell)
            assert all(math.isfinite(float(v)) for key, v in row.items() if key != 'cell')
            assert all(float(row[key]) > 0 for key in ['mass', 'volume', 'density', 'area_mean_pressure', 'U', 'K_point', 'inertia'])
    digest(f)
assert ids == cell_ids

headers = []
for rank in range(4):
    f = STEP / f'accepted_r{rank}.bin'
    with f.open('rb') as stream:
        header = struct.unpack('<8Q', stream.read(64))
    assert header[0] == 0x4d4645564f4c3031 and header[1:4] == (3, rank, 4)
    assert header[5] == 32 and header[7] == 16
    headers.append({'rank': rank, 'header': header, 'bytes': f.stat().st_size, 'sha256': digest(f)})

assert cost['projection_is_optimistic'] and cost['future_guard_cost_and_mesh_lifetime_unproved']
selected = cost['selected_step_seconds_excluding_separately_timed_reference']
h = cost['sustainable_dt_at_current_original_node_bound']
projected = (1 - result['accepted_time']) * selected / h
assert math.isclose(projected, cost['projected_remaining_seconds_without_future_guard_overhead'], rel_tol=1e-8)
speedup = 225 / selected
assert speedup < 10 and projected > 86400
assert result['stop_reason'] == 'measured_cost_exceeds_24h_t1_budget'
assert not contract['physical_long_term_heating_resolved'] and not contract['long_Kepler_validated']
summary = {
    'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'audit': 'PASS_COMPLETED_ACTUAL_PAIRED_STEP_COST_GATE_FAILED',
    'scope': 'Read completed native records, scan all low-stage and endpoint cell summaries, verify all rank/pass completion and hash accepted shards; no evolution or independent recomputation of point forces.',
    'integrity_records_all_OK': integrity,
    'same_trajectory_accepted_time': result['accepted_time'], 'total_accepted_steps': 5,
    'high_k1_force_and_solve_reused': True,
    'step_cached_initial_k1_false_means_low_rate_computed': True,
    'checkpoint_headers_and_SHA256': headers,
    'restart_scope': 'Actual native binary80/FNV reader checked complete checkpoint payload and compared state/ledgers; this retrospective audit checks headers and freezes SHA256 without replay.',
    'stages': stages, 'endpoint': endpoint, 'step': step, 'cost': cost,
    'speedup_against_225_second_baseline': speedup,
    'optimistic_t1_days': projected / 86400,
    'rank0_low_stage_mass_solve_seconds_sum': sum(item['rhs']['mass_solve_seconds_rank0'] for item in stages),
    'rank0_low_stage_prepare_seconds_sum': sum(item['rhs']['prepare_seconds_rank0'] for item in stages),
    'cost_strategy': 'Stop further fixed-10.8M-DOF optimization and blind long continuation under adopted Astra decision; preserve accepted step5 and request one strategic replacement plan.',
    't6_10_stability_proven': False, 'longest_GIZMO_stability_proven': False,
    'global_injectivity_proven': False, 'physical_long_term_heating_resolved': False,
    'no_new_calculation_launched': True,
}
(DEST / 'contract.json').write_text(json.dumps(summary, indent=2) + '\n')
digest(DEST / 'contract.json')
digest(Path(__file__).resolve())
(DEST / 'evidence.sha256').write_text(''.join(f'{sha}  {name}\n' for name, sha in sorted(evidence.items())))
print(json.dumps({key: summary[key] for key in ['audit', 'same_trajectory_accepted_time', 'optimistic_t1_days', 'speedup_against_225_second_baseline', 'rank0_low_stage_mass_solve_seconds_sum', 'rank0_low_stage_prepare_seconds_sum']}, indent=2))
