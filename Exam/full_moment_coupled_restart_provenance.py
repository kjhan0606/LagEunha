"""Seal and verify accepted self-consistent-midpoint checkpoints.

This module authenticates an already accepted state and its audits. It does
not update, retry, or accept a physical step. Each seal permits only the next
predeclared time gate at the same h.
"""
import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from full_moment_restart_provenance import (
    operator_identity,
    physics_sha,
    sha,
    verify_source_runtime_identity,
)
from full_moment_state_io import FIELDS, HEADER, digest as state_digest, dump, load

INTEGRATOR = 'self-consistent midpoint + Sylvester/exponential/positive QR'
FORMAT = 'full_moment_native_coupled_restart_provenance_v1'
STAGES = ((.03, .05),) + tuple((i/40, (i+1)/40) for i in range(2, 400))
EPS = np.finfo(float).eps


def next_gate(time):
    for source, target in STAGES:
        if math.isclose(time, source, rel_tol=0., abs_tol=32*EPS*max(1., abs(source))):
            return target
    raise AssertionError(f'No predeclared coupled continuation gate from t={time!r}')


def assert_roundoff_equal(a, b, label):
    assert math.isfinite(float(a)) and math.isfinite(float(b))
    assert abs(float(a)-float(b)) <= 32*EPS*max(1., abs(float(a)), abs(float(b))), label


def source_identity(checkpoint, expected_ranks):
    identity_path, identity, operator_sha, sources, params_path = \
        verify_source_runtime_identity(checkpoint, expected_ranks)
    provider_script = Path(identity['provider_script']).resolve(strict=True)
    provenance_path = provider_script.parent/'provenance.json'
    provenance = json.loads(provenance_path.read_text())
    assert provenance['binary_sha256'] == identity['binary_sha256']
    assert provenance['sources'] == sources
    build_manifest_path = Path(provenance['binary_build_provenance_file']).resolve(strict=True)
    assert sha(build_manifest_path) == provenance['binary_build_provenance_sha256']
    build = json.loads(build_manifest_path.read_text())
    assert build['binary_sha256'] == identity['binary_sha256']
    assert build['user_main_sha256'] == provenance['user_source_sha256']
    return (identity_path, identity, operator_sha, sources, params_path,
            provenance_path, build_manifest_path, build)


def seal(checkpoint, initial, physics, pair_gate, output):
    checkpoint, initial, physics, output = map(
        lambda p: Path(p).resolve(), (checkpoint, initial, physics, output))
    pair_gate = Path(pair_gate).resolve()
    assert checkpoint.is_dir() and checkpoint.name.startswith('accepted_')
    state = load(checkpoint)
    header_bytes = (checkpoint/'header.bin').read_bytes()
    header = HEADER.unpack(header_bytes)
    step = int(checkpoint.name.split('_', 1)[1])
    time = float(state['time'])
    assert step > 0 and time > 0 and header[3] > 0 and header[4] > 0

    parent = checkpoint.parent
    accepted_path = parent/f'accepted_{step:04d}.json'
    completed_path = parent/f'completed_{step:04d}.json'
    accepted = json.loads(accepted_path.read_text())
    completed = json.loads(completed_path.read_text())
    for report in (accepted, completed):
        assert report['model'] == 'FullMomentKepler'
        assert report['integrator'] == INTEGRATOR
        assert report['accepted_commit'] and report['fixed_mass'] and report['Re_infinity']
        assert report['step'] == step and report['cells'] == header[2]
        assert report['ranks'] in (1, 2) and report['h'] in (.0025, .005)
        assert_roundoff_equal(report['time'], time, 'Report/checkpoint time mismatch')
        assert_roundoff_equal(report['time'], step*report['h'], 'Global time-step mismatch')
    assert accepted['status'] == 'accepted_native_coupled_midpoint_step'
    assert completed['status'] == 'completed_native_evolution_batch'
    assert completed['h'] == accepted['h'] and completed['ranks'] == accepted['ranks']

    control_path = parent/f'coupled_{step:04d}'/'control.json'
    control = json.loads(control_path.read_text())
    assert control['status_code'] == 0 and control['converged'] is True
    assert control['final_checks_passed'] is True and control['outer_native_commit_required'] is True
    assert control['provider_calls'] > 0 and control['accepted'] is False

    audit = parent/f'candidate_audit_{step:04d}'
    admissibility_path = audit/'admissibility.json'
    admissibility = json.loads(admissibility_path.read_text())
    assert admissibility['status'] == 'passed_native_tentative_admissibility'
    assert admissibility['mass_bitwise_initial'] and admissibility['minimum_thermal'] > 0
    assert all(abs(x) <= bound for x, bound in zip(
        admissibility['P_E_L_ledger_residual'], admissibility['P_E_L_ledger_bound']))
    assert abs(admissibility['combined_energy_change']) <= admissibility['combined_energy_roundoff_bound']

    physical_path = audit/'physical_gate.json'
    physical = json.loads(physical_path.read_text())
    assert physical['status'] == 'passed_checkpoint_physical_monitor'
    assert physical['regions'] and all(x['passed'] for x in physical['regions'].values())

    density_path = audit/'density'/'report.json'
    density = json.loads(density_path.read_text())
    candidate_path = audit/'candidate.npz'
    assert density['status'] == 'passed_full_stage_Psi_density_preparation'
    assert density['source'] == 'Native-owned fits with unchanged independent q48'
    assert density['cells'] == header[2] and density['time'] == time
    assert sha(candidate_path) == density['state_sha256']
    candidate = load(candidate_path)
    assert all(np.array_equal(state[k], candidate[k]) for k in state)

    agreement_path = parent/f'coupled_{step:04d}'/'agreement'/'report.json'
    agreement = json.loads(agreement_path.read_text())
    assert agreement['status'] == 'passed_actual_midpoint_agreement'
    journal_path = parent/f'coupled_{step:04d}'/'accepted_journal'/'report.json'
    journal = json.loads(journal_path.read_text())
    assert journal['status'] == 'captured_verified_map_cell_journal'
    assert_roundoff_equal(journal['time'], time, 'Coupled journal/checkpoint time mismatch')
    assert_roundoff_equal(journal['h'], accepted['h'], 'Coupled journal/checkpoint h mismatch')

    old = load(initial)
    assert float(old['time']) == 0.
    assert not state['frozen'].any() and state['psi_constrained_mask'].all()
    assert np.array_equal(state['ids'], old['ids'])
    assert np.array_equal(state['mass'], old['mass'])

    pair = json.loads(pair_gate.read_text())
    assert pair['status'] == 'passed_coupled_pair_gate'
    assert pair['all_physical_gates_passed'] is True
    assert pair['fine_coarse_consistency_passed'] is True
    assert_roundoff_equal(pair['time'], time, 'Pair gate/checkpoint time mismatch')
    assert set(pair['branches']) == {'fine', 'coarse'}
    comparison_path = Path(pair['comparison_artifact_path']).resolve(strict=True)
    analysis_script_path = Path(pair['analysis_script_path']).resolve(strict=True)
    assert sha(comparison_path) == pair['comparison_artifact_sha256']
    assert sha(analysis_script_path) == pair['analysis_script_sha256']
    matched_branch = None
    for label, branch in pair['branches'].items():
        branch_path = Path(branch['checkpoint_path']).resolve()
        assert branch['accepted'] and branch['physical_gates_passed']
        assert branch['integrator'] == INTEGRATOR
        assert branch['h'] == (.0025 if label == 'fine' else .005)
        assert branch['step'] > 0
        assert_roundoff_equal(branch['step']*branch['h'], time, 'Pair branch global step mismatch')
        assert branch['checkpoint_sha256'] == state_digest(branch_path)
        assert_roundoff_equal(branch['time'], time, 'Fine/coarse pair endpoint mismatch')
        if branch_path == checkpoint:
            matched_branch = label
    assert matched_branch is not None, 'This exact accepted checkpoint was not included in the pair gate'

    runtime_identity_path, identity, source_operator_sha, provider_sources, params_path, \
        provider_provenance_path, build_manifest_path, build = source_identity(checkpoint, accepted['ranks'])
    operator, operator_sha = operator_identity()
    assert source_operator_sha == operator_sha and build['user_main_sha256']
    for rank in range(accepted['ranks']):
        factor_path = parent/f'factor_{header[3]}_rank{rank}.json'
        factor = json.loads(factor_path.read_text())
        assert factor['status'] == 'passed_owned_native_factor_audit'
        assert factor['rank'] == rank and factor['minimum_thermal'] > 0
        assert factor['maximum_QR_bound_fraction'] <= 1 and factor['shadow_audit']

    limit = next_gate(time)
    evidence_paths = [accepted_path, completed_path, control_path, admissibility_path,
                      physical_path, density_path, candidate_path, agreement_path,
                      journal_path, pair_gate, comparison_path, analysis_script_path,
                      runtime_identity_path, params_path,
                      provider_provenance_path, build_manifest_path]
    evidence_paths.extend(parent/f'factor_{header[3]}_rank{rank}.json'
                          for rank in range(accepted['ranks']))
    state_files = {name+'.bin': sha(checkpoint/(name+'.bin')) for name in ('header', *FIELDS)}
    record = dict(format=FORMAT, model='FullMomentKepler', integrator=INTEGRATOR,
        fixed_mass=True, Re_infinity=True, steering=False,
        source=str(checkpoint), source_sha256=state_digest(checkpoint), time=time, step=step,
        h=accepted['h'], ranks=accepted['ranks'], cells=header[2],
        original_epoch=header[3], original_run_token=header[4],
        next_gate_time=limit, maximum_target_time=limit, allowed_target_h=[accepted['h']],
        pair_gate_path=str(pair_gate), pair_gate_sha256=sha(pair_gate),
        initial_state_sha256=sha(initial), physics_sha256=physics_sha(physics),
        operator_path=str(operator), operator_sha256=operator_sha,
        source_runtime_identity_path=str(runtime_identity_path),
        source_runtime_identity_sha256=sha(runtime_identity_path),
        source_params_path=str(params_path), source_params_sha256=sha(params_path),
        source_binary_sha256=identity['binary_sha256'],
        source_user_main_sha256=build['user_main_sha256'],
        source_operator_sha256=source_operator_sha,
        source_provider_script=identity['provider_script'],
        source_provider_snapshot_sources_sha256=provider_sources,
        source_provider_provenance_sha256=sha(provider_provenance_path),
        source_build_manifest_sha256=sha(build_manifest_path),
        state_files=state_files,
        evidence={str(path): sha(path) for path in evidence_paths},
        physical_providers=0, accepted_time_steps=0)
    output.mkdir(exist_ok=False)
    path = output/'manifest.json'
    with path.open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    return path


def verify_serialize(source, manifest, initial, physics, requested_h,
                     target_step, maximum_target_time, output):
    source, manifest, initial, physics, output = map(
        lambda p: Path(p).resolve(), (source, manifest, initial, physics, output))
    record = json.loads(manifest.read_text())
    assert record['format'] == FORMAT and record['source'] == str(source)
    assert record['source_sha256'] == state_digest(source)
    assert record['model'] == 'FullMomentKepler' and record['integrator'] == INTEGRATOR
    assert record['fixed_mass'] and record['Re_infinity'] and record['steering'] is False
    assert record['initial_state_sha256'] == sha(initial)
    assert record['physics_sha256'] == physics_sha(physics)
    assert math.isfinite(requested_h) and requested_h == record['h']
    assert requested_h in record['allowed_target_h']
    assert math.isfinite(maximum_target_time)
    assert maximum_target_time == record['maximum_target_time'] == record['next_gate_time']
    target_time = target_step*requested_h
    assert target_step > record['step']
    assert target_time > record['time'] and target_time <= maximum_target_time
    assert_roundoff_equal(target_time, target_step*record['h'], 'Requested target time mismatch')

    operator, operator_sha = operator_identity()
    assert record['operator_path'] == str(operator) and record['operator_sha256'] == sha(operator) == operator_sha
    identity_path, identity, source_operator_sha, sources, params_path, \
        provider_provenance_path, build_manifest_path, build = source_identity(source, record['ranks'])
    assert record['source_runtime_identity_path'] == str(identity_path)
    assert record['source_runtime_identity_sha256'] == sha(identity_path)
    assert record['source_params_path'] == str(params_path)
    assert record['source_params_sha256'] == sha(params_path) == identity['params_sha256']
    assert record['source_binary_sha256'] == identity['binary_sha256']
    assert record['source_user_main_sha256'] == build['user_main_sha256']
    assert record['source_operator_sha256'] == source_operator_sha
    assert record['source_provider_script'] == identity['provider_script']
    assert record['source_provider_snapshot_sources_sha256'] == sources
    assert record['source_provider_provenance_sha256'] == sha(provider_provenance_path)
    assert record['source_build_manifest_sha256'] == sha(build_manifest_path)

    for name, digest in record['state_files'].items():
        assert sha(source/name) == digest, f'Checkpoint changed: {name}'
    for path, digest in record['evidence'].items():
        assert sha(path) == digest, f'Coupled restart evidence changed: {path}'
    header = (source/'header.bin').read_bytes()
    fields = HEADER.unpack(header)
    assert fields[2] == record['cells'] and fields[3] == record['original_epoch']
    assert fields[4] == record['original_run_token'] and fields[8] == record['time']
    assert fields[8] == record['step']*record['h']
    state = load(source)
    old = load(initial)
    assert np.array_equal(state['ids'], old['ids']) and np.array_equal(state['mass'], old['mass'])

    output.mkdir(exist_ok=False)
    token = int(record['source_sha256'][:16], 16) or 1
    serialized = dump(output/'serialized_source', state, fields[3], token)
    expected_header = (serialized/'header.bin').read_bytes()
    with (output/'expected_header.bin').open('xb') as stream:
        stream.write(expected_header)
        stream.flush()
        os.fsync(stream.fileno())
    observed = load(serialized)
    assert all(np.array_equal(observed[k], state[k]) for k in state)
    receipt = dict(status='verified_coupled_restart_provenance_no_state_update',
        source=str(source), source_sha256=record['source_sha256'],
        manifest_sha256=sha(manifest), source_integrator=INTEGRATOR,
        target_integrator=INTEGRATOR, source_time=record['time'], source_step=record['step'],
        target_time=target_time, target_step=target_step, maximum_target_time=maximum_target_time,
        h=requested_h, cells=record['cells'], physics_sha256=record['physics_sha256'],
        serialized_source=str(serialized), serialization_token=token,
        payload_sha256={name+'.bin': sha(serialized/(name+'.bin')) for name in FIELDS},
        all_six_payload_values_preserved=True, physical_providers=0, accepted_time_steps=0)
    with (output/'receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    return expected_header, receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--initial', type=Path, required=True)
    parser.add_argument('--physics', type=Path, required=True)
    parser.add_argument('--pair-gate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(seal(args.checkpoint, args.initial, args.physics, args.pair_gate, args.output), flush=True)


if __name__ == '__main__':
    main()
