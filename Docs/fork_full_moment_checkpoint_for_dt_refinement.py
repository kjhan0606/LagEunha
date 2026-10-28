#!/usr/bin/env python3
"""Create explicit, provenance-bound timestep-fork inputs from an accepted NPZ.

This is data serialization and audit plumbing only: it performs no RHS,
geometry update, factor update, or timestep. The target h is explicit in each
fork manifest and remains guarded by the native restart verifier.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

BASE = Path('/gpfs/kjhan/LagForce')
SOURCE_ROOT = BASE / 'Codex_full_disk_continuous_two_dt_to_t1_20261019'
SOURCE_FINE = SOURCE_ROOT / 'fine'
SOURCE_CHECKPOINT = SOURCE_FINE / 'checkpoints/fine_0010'
CHECKPOINT_NPZ = SOURCE_CHECKPOINT / 'accepted.npz'
INITIAL = BASE / 'Codex_full_curved_kepler_IC_20261019/initial.npz'
PHYSICS = BASE / 'Codex_full_curved_kepler_IC_20261019/report.json'
BUILD = BASE / 'Codex_full_disk_native_evolution_build_v10_20261019'
EXE = BUILD / 'eunha2.exe'
PROVIDER = BUILD / 'full_moment_provider_service.py'
RUNNER = BASE / 'Codex_full_disk_native_restart_equivalence_v10_20261019/run_full_moment_with_identity.py'
OPERATOR = BASE / 'Codex_full_disk_owned_ledger_build_v4_20261019/liblagmoment_owned.so'
OUTPUT = BASE / 'Codex_full_disk_explicit_fine_timestep_refinement_fork_v2_20261019'
EXPECTED_CHECKPOINT_SHA = '465f9713ab405f19e438396abb8796510c3c51784b0c85028139248b34c28d8e'
EXPECTED_INITIAL_SHA = '6794d8e56c0bf8829bf95bb2ae1a633d2b690588db243ccbccdf9f8d45809c46'
EXPECTED_OPERATOR_SHA = '51014e2fb7de34bb8eb1ed99c8ed0e22ecbbe5690a67327f8fec9892d5968471'
EXPECTED_PHYSICS_SHA = '5128e9d2475f60389ba6b59aaeab34a64d5f75d96bbaf91fa9d13ee9f8d55efc'
SOURCE_TIME = 0.025
SOURCE_H = 0.0025
SOURCE_STEP = 10
END_TIME = 0.05
TARGETS = ((0.0025, 'h0025', '06:00:00'), (0.00125, 'h00125', '10:00:00'))

sys.path.insert(0, str(BUILD))
from full_moment_state_io import FIELDS, dump, load  # noqa: E402
from full_moment_restart_provenance import (  # noqa: E402
    FORMAT,
    INTEGRATOR,
    OPERATOR_SHA,
    physics_sha,
    verify,
)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block := stream.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def write_json_exclusive(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def check_source() -> tuple[dict, dict, dict, dict, list[Path]]:
    assert not OUTPUT.exists(), f'refusing to overwrite {OUTPUT}'
    assert sha(CHECKPOINT_NPZ) == EXPECTED_CHECKPOINT_SHA
    assert sha(INITIAL) == EXPECTED_INITIAL_SHA
    assert sha(OPERATOR) == EXPECTED_OPERATOR_SHA == OPERATOR_SHA
    assert sha(EXE) == (BUILD / 'binary.sha256').read_text().split()[0]
    assert RUNNER.is_file()

    fine_manifest = read_json(SOURCE_FINE / 'manifest.json')
    fine_result = read_json(SOURCE_FINE / 'result.json')
    assert fine_manifest['status'] == 'starting_approved_continuous_bundle'
    assert fine_manifest['mode'] == 'fine'
    assert fine_manifest['fixed_mass'] and fine_manifest['Re'] == 'infinity'
    assert fine_manifest['physical_viscosity'] == 0 and fine_manifest['steering'] is False
    assert fine_manifest['initial_sha256'] == sha(INITIAL)
    assert fine_manifest['cells'] == 19945
    assert fine_manifest['steps']['fine_h'] == SOURCE_H
    assert fine_manifest['steps']['comparison_end'] == END_TIME
    assert fine_result['status'] == 'completed_comparison_branch'
    assert fine_result['branch'] == 'fine' and fine_result['time'] == END_TIME
    assert fine_result['accepted_steps'] == 20

    state = load(CHECKPOINT_NPZ)
    assert float(state['time']) == SOURCE_TIME and len(state['ids']) == 19945
    initial = load(INITIAL)
    assert np.array_equal(state['ids'], initial['ids'])
    assert np.array_equal(state['mass'], initial['mass'])
    assert not state['frozen'].any() and state['psi_constrained_mask'].all()

    gate = read_json(SOURCE_CHECKPOINT / 'physical_gate.json')
    state_checks = read_json(SOURCE_CHECKPOINT / 'state_checks.json')
    density_prepare = read_json(SOURCE_CHECKPOINT / 'density_prepare.json')
    density = read_json(SOURCE_CHECKPOINT / 'density/report.json')
    physics = read_json(SOURCE_CHECKPOINT / 'physics/report.json')
    reference = read_json(SOURCE_CHECKPOINT / 'reference/report.json')
    assert gate['status'] == 'passed_checkpoint_physical_monitor'
    assert gate['time'] == SOURCE_TIME and all(x['passed'] for x in gate['regions'].values())
    assert state_checks['mass_bitwise_initial'] and state_checks['minimum_thermal'] > 0
    assert all(abs(a) <= b for a, b in zip(
        state_checks['P_E_L_ledger_residual'], state_checks['P_E_L_ledger_bound']))
    assert abs(state_checks['combined_energy_change']) <= state_checks['combined_energy_roundoff_bound']
    assert density_prepare['status'] == 'passed_in_memory_Psi_density_preparation'
    assert density_prepare['independent_q48'] and density_prepare['cells'] == len(state['ids'])
    assert density['status'] == 'passed_full_stage_Psi_density_preparation'
    assert density['time'] == SOURCE_TIME and density['state_sha256'] == EXPECTED_CHECKPOINT_SHA
    assert physics['evolution_steps'] == 0 and physics['provider_calls'] == 0
    assert reference['actual_state_sha256'] == EXPECTED_CHECKPOINT_SHA
    assert reference['counts']['passed'] == len(state['ids'])
    assert read_json(PHYSICS)['physics']['Re'] == 'infinity'
    assert physics_sha(PHYSICS) == EXPECTED_PHYSICS_SHA

    provenance = read_json(BUILD / 'provenance.json')
    assert provenance['binary_sha256'] == sha(EXE)
    for name, digest in provenance['sources'].items():
        assert sha(BUILD / name) == digest, f'build source changed: {name}'

    evidence = [
        SOURCE_FINE / 'manifest.json', SOURCE_FINE / 'result.json',
        CHECKPOINT_NPZ, SOURCE_CHECKPOINT / 'physical_gate.json',
        SOURCE_CHECKPOINT / 'state_checks.json', SOURCE_CHECKPOINT / 'density_prepare.json',
        SOURCE_CHECKPOINT / 'density/report.json', SOURCE_CHECKPOINT / 'physics/report.json',
        SOURCE_CHECKPOINT / 'reference/report.json', INITIAL, PHYSICS,
        BUILD / 'provenance.json', BUILD / 'binary.sha256', RUNNER,
    ]
    return state, gate, density, reference, evidence


def submit_script(run: Path, params: Path, target_h: float, walltime: str, label: str) -> str:
    return f'''#!/bin/bash
#SBATCH --job-name=kep-{label}
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=128G
#SBATCH --time={walltime}
#SBATCH --chdir={run}
#SBATCH --output={run}/slurm-%j.out
set -euo pipefail
module purge
module load intel/tbb intel/compiler-rt intel/umf intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2
export LD_LIBRARY_PATH=/home/kjhan/local/lib:/opt/ohpc/pub/compiler/gcc/13.2.0/lib64:${{LD_LIBRARY_PATH:-}}
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 I_MPI_PIN_DOMAIN=omp
export EUNHA_MOMENT_BINARY_SHA={sha(EXE)}
/home/kjhan/miniconda3/bin/python3 {RUNNER} \\
  --exe {EXE} --provider {PROVIDER} --params {params} --ranks 1 \\
  --log {run}/runtime/native.log --identity {run}/runtime_identity.json \\
  > {run}/launch.log 2>&1
'''


def main() -> None:
    state, gate, density, reference, evidence = check_source()
    OUTPUT.mkdir(exist_ok=False)
    builder_copy = OUTPUT / 'fork_full_moment_checkpoint_for_dt_refinement.py'
    with builder_copy.open('xb') as stream:
        stream.write(Path(__file__).read_bytes())
        stream.flush()
        os.fsync(stream.fileno())

    token = int.from_bytes(hashlib.sha256(
        (EXPECTED_CHECKPOINT_SHA + 'explicit-fine-timestep-fork').encode()).digest()[:8], 'little') or 1
    native_source = OUTPUT / 'source_state'
    dump(native_source, state, epoch=SOURCE_STEP, token=token)
    source_state_digest = {
        name: sha(native_source / name)
        for name in ('header.bin', *(key + '.bin' for key in FIELDS))
    }
    evidence_hashes = {str(path.resolve()): sha(path) for path in evidence}

    targets = []
    for target_h, label, walltime in TARGETS:
        start_step = round(SOURCE_TIME / target_h)
        end_step = round(END_TIME / target_h)
        assert start_step * target_h == SOURCE_TIME
        assert end_step * target_h == END_TIME
        run = OUTPUT / f'run_{label}'
        manifest_dir = run / 'input_manifest'
        runtime = run / 'runtime'
        manifest_dir.mkdir(parents=True, exist_ok=False)
        runtime.mkdir(exist_ok=False)

        manifest = dict(
            format=FORMAT, model='FullMomentKepler', integrator=INTEGRATOR,
            fixed_mass=True, Re_infinity=True, steering=False,
            source=str(native_source.resolve()), time=SOURCE_TIME,
            step=start_step, h=target_h, cells=len(state['ids']),
            original_epoch=SOURCE_STEP, original_run_token=token,
            initial_state_sha256=sha(INITIAL), physics_sha256=physics_sha(PHYSICS),
            operator_sha256=sha(OPERATOR), state_files=source_state_digest,
            evidence=evidence_hashes, physical_providers=0, accepted_time_steps=0,
            fork=dict(
                status='explicitly_authorized_timestep_refinement_branch',
                source_checkpoint=str(CHECKPOINT_NPZ.resolve()),
                source_checkpoint_sha256=EXPECTED_CHECKPOINT_SHA,
                source_integrator='explicit midpoint + Sylvester/exponential/positive QR',
                source_h=SOURCE_H, source_step=SOURCE_STEP, source_time=SOURCE_TIME,
                target_h=target_h, target_start_grid_index=start_step,
                target_end_time=END_TIME, target_end_grid_index=end_step,
                exact_payload_copy=True, no_rhs_or_state_update_during_fork=True,
            ),
        )
        manifest_path = manifest_dir / 'manifest.json'
        write_json_exclusive(manifest_path, manifest)
        header, checked = verify(native_source, manifest_path, INITIAL, PHYSICS, target_h)
        assert header and checked['fork']['target_h'] == target_h

        params = run / 'params.dat'
        params_text = f'''define Simulation Model = FullMomentKepler
moment_native_evolution = 1
moment_python = /home/kjhan/miniconda3/bin/python3
moment_provider_script = {PROVIDER}
moment_state = {native_source}
moment_initial = {INITIAL}
moment_geometry = {BASE}/Codex_full_disk_native_initial_geometry_20261019
moment_output = {runtime}
moment_physics_file = {PHYSICS}
moment_physics_sha = {EXPECTED_PHYSICS_SHA}
moment_workers = 7
moment_threads = 8
moment_h = {target_h:.8g}
moment_steps = {end_step}
moment_validation = 0
moment_restart_manifest = {manifest_path}
'''
        with params.open('x') as stream:
            stream.write(params_text)
            stream.flush()
            os.fsync(stream.fileno())
        submit = run / 'submit.sh'
        with submit.open('x') as stream:
            stream.write(submit_script(run, params, target_h, walltime, label))
            stream.flush()
            os.fsync(stream.fileno())
        submit.chmod(0o750)

        targets.append(dict(label=label, run=str(run), h=target_h,
            start_step=start_step, end_step=end_step, walltime=walltime,
            gpu_request='--gres=gpu:H200:1', ranks=1,
            manifest_sha256=sha(manifest_path), params_sha256=sha(params),
            submit_sha256=sha(submit)))

    write_json_exclusive(OUTPUT / 'fork_bundle.json', dict(
        status='passed_source_audit_and_native_restart_manifest_verification',
        source_checkpoint=str(CHECKPOINT_NPZ.resolve()),
        source_checkpoint_sha256=EXPECTED_CHECKPOINT_SHA,
        source_time=SOURCE_TIME, source_h=SOURCE_H,
        source_physical_gate=gate['status'], source_density_gate=density['status'],
        source_actual_cell_reference_passed=reference['counts']['passed'],
        native_source=str(native_source.resolve()), native_source_payload_sha256=source_state_digest,
        native_source_token=token, original_source_unchanged=True,
        physical_providers=0, accepted_time_steps=0,
        build_binary_sha256=sha(EXE), operator_sha256=sha(OPERATOR),
        runtime_identity_runner=str(RUNNER), runtime_identity_runner_sha256=sha(RUNNER),
        builder_sha256=sha(builder_copy), targets=targets,
        comparison_note='Both branches begin from the exact same fine t=.025 accepted NPZ state. The h=.0025 branch is an H200-native control; h=.00125 is the third timestep. Neither branch replays t=0 to .025.',
    ))
    print(json.dumps(dict(status='prepared_dt_refinement_fork_inputs',
        output=str(OUTPUT), targets=targets, native_source=str(native_source),
        physical_providers=0, accepted_time_steps=0), indent=2), flush=True)


if __name__ == '__main__':
    main()
