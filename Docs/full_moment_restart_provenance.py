"""Bind a completed native checkpoint to its model, baseline and audit files.

This is accidental-corruption/provenance checking, not a security signature.
No state update, geometry replacement, retry, or physical acceptance is done.
Existing v7/v8 checkpoints are sealed in a NEW directory, never rewritten.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
from full_moment_state_io import HEADER, FIELDS, load

BASE = Path('/gpfs/kjhan/LagForce')
DEFAULT_OPERATOR = BASE/'Codex_full_disk_owned_ledger_build_v4_20261019/liblagmoment_owned.so'
DEFAULT_OPERATOR_SHA = '51014e2fb7de34bb8eb1ed99c8ed0e22ecbbe5690a67327f8fec9892d5968471'
INTEGRATOR = 'explicit midpoint + Sylvester/exponential/positive QR'
FORMAT = 'full_moment_native_restart_provenance_v1'


def sha(file):
    h = hashlib.sha256()
    with Path(file).open('rb') as f:
        while block := f.read(1024*1024):
            h.update(block)
    return h.hexdigest()


def physics_sha(file):
    physics = json.loads(Path(file).read_text())['physics']
    return hashlib.sha256(json.dumps(physics, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def operator_identity():
    """Resolve an explicitly selected, content-pinned native CUDA operator."""
    operator = Path(os.environ.get('EUNHA_MOMENT_OPERATOR_PATH', str(DEFAULT_OPERATOR))).resolve()
    expected = os.environ.get('EUNHA_MOMENT_OPERATOR_SHA', DEFAULT_OPERATOR_SHA)
    actual = sha(operator)
    assert actual == expected, f'Operator SHA mismatch for {operator}: {actual} != {expected}'
    return operator, actual


def verify_source_runtime_identity(checkpoint, expected_ranks):
    """Bind a new seal to the completed run that produced its checkpoint."""
    path = checkpoint.parent.parent/'runtime_identity.json'
    identity = json.loads(path.read_text())
    assert identity['status'] == 'passed_live_native_runtime_identity'
    assert identity['native_exit_code'] == 0 and not identity['identity_errors']
    assert identity['binary_sha256'] and len(identity['binary_sha256']) == 64
    ranks = {str(rank) for rank in range(expected_ranks)}
    assert set(identity['ranks']) == ranks
    assert set(identity['GPU_UUID_by_rank']) == ranks
    assert len(set(identity['GPU_UUID_by_rank'].values())) == expected_ranks

    executable_paths = set()
    operator_hashes = set()
    provider_script = Path(identity['provider_script']).resolve(strict=True)
    for rank in sorted(ranks, key=int):
        row = identity['ranks'][rank]
        executable = Path(row['actual_executable']).resolve(strict=True)
        assert sha(executable) == identity['binary_sha256']
        executable_paths.add(str(executable))
        libraries = row['actual_mapped_operator_libraries']
        assert len(libraries) == 1
        for library, digest in libraries.items():
            library = Path(library).resolve(strict=True)
            assert sha(library) == digest
            operator_hashes.add(digest)
        providers = row['actual_provider_process']
        assert len(providers) == 1
        assert any(str(provider_script) in argv for argv in providers.values())
    assert len(executable_paths) == 1 and len(operator_hashes) == 1

    sources = identity['provider_snapshot_sources_sha256']
    assert sources
    for name, digest in sources.items():
        assert Path(name).name == name, f'Invalid provider source path: {name}'
        assert sha(provider_script.parent/name) == digest, f'Provider source changed: {name}'
    params = path.parent/'params.dat'
    assert sha(params) == identity['params_sha256'], 'Run parameters changed since launch'
    return path, identity, next(iter(operator_hashes)), sources, params


def seal(checkpoint, initial, physics, output):
    """Require existing native acceptance, checkpoint q48 and physical audits."""
    checkpoint, initial, physics, output = map(lambda p: Path(p).resolve(),
                                              (checkpoint, initial, physics, output))
    state = load(checkpoint)
    header_bytes = (checkpoint/'header.bin').read_bytes()
    header = HEADER.unpack(header_bytes)
    assert header[3] and header[4] and float(state['time']) > 0
    operator, operator_sha = operator_identity()
    assert checkpoint.name.startswith('accepted_')
    step = int(checkpoint.name.split('_')[1])
    parent = checkpoint.parent
    journal = parent/f'accepted_{step:04d}.json'
    report = json.loads(journal.read_text())
    assert report['status'] == 'accepted_native_midpoint_step'
    assert report['native_time_control'] and report['accepted_commit']
    assert report['model'] == 'FullMomentKepler' and report['integrator'] == INTEGRATOR
    assert report['fixed_mass'] and report['Re_infinity']
    assert report['step'] == step and report['cells'] == header[2]
    assert report['time'] == header[8] == step*report['h']
    assert report['ranks'] in (1, 2)
    runtime_identity_path, runtime_identity, source_operator_sha, source_provider_sources, runtime_params_path = \
        verify_source_runtime_identity(checkpoint, report['ranks'])
    old = load(initial)
    assert np.array_equal(state['ids'], old['ids']) and np.array_equal(state['mass'], old['mass'])
    assert float(old['time']) == 0 and not state['frozen'].any() and state['psi_constrained_mask'].all()
    audit = parent/f'candidate_audit_{step:04d}'
    admissibility = json.loads((audit/'admissibility.json').read_text())
    assert admissibility['status'] == 'passed_native_tentative_admissibility'
    assert admissibility['time'] == header[8] and admissibility['epoch'] == header[3]
    assert admissibility['mass_bitwise_initial'] and admissibility['minimum_thermal'] > 0
    assert all(abs(a) <= b for a, b in zip(admissibility['P_E_L_ledger_residual'], admissibility['P_E_L_ledger_bound']))
    assert abs(admissibility['combined_energy_change']) <= admissibility['combined_energy_roundoff_bound']
    physical = json.loads((audit/'physical_gate.json').read_text())
    assert physical['status'] == 'passed_checkpoint_physical_monitor' and physical['time'] == header[8]
    assert physical['regions'] and all(row['passed'] for row in physical['regions'].values())
    density = json.loads((audit/'density/report.json').read_text())
    assert density['status'] == 'passed_full_stage_Psi_density_preparation' and density['cells'] == header[2]
    assert density['time'] == header[8] and density['source'] == 'Native-owned fits with unchanged independent q48'
    candidate = audit/'candidate.npz'
    assert sha(candidate) == density['state_sha256']
    candidate_state = load(candidate)
    assert all(np.array_equal(state[k], candidate_state[k]) for k in state)
    evidence = [journal, audit/'admissibility.json', audit/'physical_gate.json',
                audit/'density/report.json', candidate, runtime_identity_path, runtime_params_path]
    for rank in range(report['ranks']):
        factor = parent/f'factor_{header[3]}_rank{rank}.json'
        row = json.loads(factor.read_text())
        assert row['status'] == 'passed_owned_native_factor_audit' and row['rank'] == rank
        assert row['minimum_thermal'] > 0 and row['maximum_QR_bound_fraction'] <= 1
        assert row['shadow_audit'], 'Restart checkpoint must have independent factor shadows'
        evidence.append(factor)
    record = dict(format=FORMAT, model='FullMomentKepler', integrator=INTEGRATOR,
        fixed_mass=True, Re_infinity=True, steering=False,
        source=str(checkpoint), time=header[8], step=step, h=report['h'], ranks=report['ranks'],
        cells=header[2], original_epoch=header[3], original_run_token=header[4],
        initial_state_sha256=sha(initial), physics_sha256=physics_sha(physics),
        operator_path=str(operator), operator_sha256=operator_sha,
        source_runtime_identity_path=str(runtime_identity_path),
        source_runtime_identity_sha256=sha(runtime_identity_path),
        source_params_path=str(runtime_params_path), source_params_sha256=sha(runtime_params_path),
        source_binary_sha256=runtime_identity['binary_sha256'],
        source_operator_sha256=source_operator_sha,
        source_provider_script=runtime_identity['provider_script'],
        source_provider_snapshot_sources_sha256=source_provider_sources,
        state_files={k+'.bin': sha(checkpoint/(k+'.bin')) for k in ('header', *FIELDS)},
        evidence={str(p): sha(p) for p in evidence},
        physical_providers=0, accepted_time_steps=0)
    output.mkdir(exist_ok=False)
    file = output/'manifest.json'
    with file.open('x') as f:
        json.dump(record, f, indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    return file


def verify(checkpoint, manifest, initial, physics, requested_h):
    checkpoint = Path(checkpoint).resolve()
    record = json.loads(Path(manifest).read_text())
    assert record['format'] == FORMAT and record['source'] == str(checkpoint)
    assert record['model'] == 'FullMomentKepler' and record['integrator'] == INTEGRATOR
    assert record['fixed_mass'] and record['Re_infinity'] and record['steering'] is False
    assert record['initial_state_sha256'] == sha(initial)
    assert record['physics_sha256'] == physics_sha(physics)
    operator, operator_sha = operator_identity()
    assert record['operator_sha256'] == sha(operator) == operator_sha
    if 'operator_path' in record:
        assert record['operator_path'] == str(operator)
    else:
        # Backward compatibility for manifests made before path pinning. Such
        # manifests can only verify against the original default operator.
        assert operator_sha == DEFAULT_OPERATOR_SHA
    if 'source_runtime_identity_path' in record:
        identity_path, identity, source_operator_sha, sources, params_path = \
            verify_source_runtime_identity(checkpoint, record['ranks'])
        assert record['source_runtime_identity_path'] == str(identity_path)
        assert record['source_runtime_identity_sha256'] == sha(identity_path)
        assert record['source_params_path'] == str(params_path)
        assert record['source_params_sha256'] == sha(params_path) == identity['params_sha256']
        assert record['source_binary_sha256'] == identity['binary_sha256']
        assert record['source_operator_sha256'] == source_operator_sha
        assert record['source_provider_script'] == identity['provider_script']
        assert record['source_provider_snapshot_sources_sha256'] == sources
    assert math.isfinite(requested_h) and requested_h == record['h'], 'No silent restart timestep change'
    assert set(record['state_files']) == {k+'.bin' for k in ('header', *FIELDS)}
    for name, digest in record['state_files'].items():
        assert sha(checkpoint/name) == digest, f'Checkpoint changed: {name}'
    assert record['evidence']
    for path, digest in record['evidence'].items():
        assert sha(path) == digest, f'Audit changed: {path}'
    header = (checkpoint/'header.bin').read_bytes()
    fields = HEADER.unpack(header)
    assert fields[2] == record['cells'] and fields[3] == record['original_epoch']
    assert fields[4] == record['original_run_token'] and fields[8] == record['time']
    assert fields[8] == record['step']*requested_h
    load(checkpoint)  # Full finite/type/count/long-double validation, no update.
    return header, record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint', type=Path, required=True)
    p.add_argument('--initial', type=Path, required=True)
    p.add_argument('--physics', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    print(seal(args.checkpoint, args.initial, args.physics, args.output), flush=True)


if __name__ == '__main__':
    main()
