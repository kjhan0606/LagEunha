"""Independent original global-stage path on the native run's H200 host.

No MPI ownership bridge. All native checkpoint fields, including the nonlinear
root warm-start seed, must match exactly. The earlier cross-host failure is
preserved; no tolerance is enlarged. CPU preparations run in fresh processes
before each process's first CUDA call. Exactly two RHS evaluations, one step.
"""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import numpy as np

BASE = Path('/gpfs/kjhan/LagForce')
FAST = BASE/'Codex_full_disk_fast_modules_v1_20261019'
OLD = BASE/'Codex_full_disk_arc_roots_modules_v7_20261019'
NATIVE = BASE/'Codex_full_disk_native_evolution_build_v8_20261019'
LIBRARY = BASE/'Codex_full_disk_owned_ledger_build_v4_20261019/liblagmoment_owned.so'
IC = BASE/'Codex_full_curved_kepler_IC_20261019/initial.npz'
sys.path.insert(0, str(FAST))
import full_disk_stage_memory as stage_api
import full_disk_density_memory as density_api
import full_disk_velocity_memory as velocity_api
from full_disk_factor_memory import update
from moment_compiled_backend import Backend
sys.path.insert(0, str(NATIVE))
from full_moment_state_io import load


def sha(file):
    return hashlib.sha256(Path(file).read_bytes()).hexdigest()


def write(file, row):
    with Path(file).open('x') as f:
        json.dump(row, f, indent=2)


def stage(args):
    state = load(args.state)
    args.output.mkdir(exist_ok=False)
    if float(state['time']) == 0:
        geometry = BASE/'Codex_full_disk_native_initial_geometry_20261019'
    else:
        geometry = args.output/'geometry'
        with (args.output/'geometry.log').open('x') as out:
            subprocess.run([sys.executable, str(OLD/'prepare_full_curved_stage_geometry_workers.py'),
                '--state', str(args.state), '--workers', '8', '--output', str(geometry)],
                check=True, stdout=out, stderr=subprocess.STDOUT)
    mesh = load(geometry/'mesh.npz')
    native = json.loads((geometry/'full_exact_supports.json').read_text())
    key = sha(args.state)+':'+sha(geometry/'full_exact_supports.json')
    fits = density_api.prepare(state, native, mesh, key, workers=8,
                               independent_q48=float(state['time']) == 0)
    write(args.output/'density.json', fits.report)
    velocity = velocity_api.prepare(state, native, mesh, fits.rows, key, workers=8)
    write(args.output/'velocity.json', velocity.report)
    backend = Backend(LIBRARY, 'cuda', 8)
    result = stage_api.evaluate(state, native, mesh, fits.rows, velocity.prepared,
        velocity.gravity_nodes, key, backend, args.output/'actual_stage', provider_budget=2,
        callback=lambda row: print(json.dumps(row), flush=True))
    np.savez(args.output/'rates.npz', **result.saved)
    write(args.report, result.report)
    del result, velocity
    gc.collect()


def main(args):
    assert socket.gethostname().split('.')[0] == 'syn104'
    args.output.mkdir(exist_ok=False)
    memory = Path(f'/dev/shm/eunha_original_oracle_{os.getpid()}_{time.time_ns()}')
    memory.mkdir(exist_ok=False)
    write(args.output/'preflight.json', dict(host=socket.gethostname(), memory=str(memory),
        original_global_stage=True, native_ownership_bridge=False, physical_provider_budget=2,
        gpu_library_sha256=sha(LIBRARY), initial_sha256=sha(IC), output_policy=
        'Two compact states plus JSON audits. Incoming factors temporary in node-local memory; retain on failure.'))
    initial = load(IC)
    backend = Backend(LIBRARY, 'cuda', 8)

    def provider(state_file, label):
        folder = memory/label
        with (args.output/(label+'.log')).open('x') as log:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), '--mode', 'stage',
                '--state', str(state_file), '--output', str(folder),
                '--report', str(args.output/(label+'.json'))], check=True,
                stdout=log, stderr=subprocess.STDOUT)
        return load(folder/'rates.npz')

    first = provider(IC, 'initial_rhs')
    midpoint, audit = update(initial, initial, first, .0005, .0005, backend, shadow_audit=True)
    write(args.output/'half_factor.json', audit)
    midpoint_file = args.output/'midpoint.npz'
    np.savez(midpoint_file, **midpoint)
    del first
    gc.collect()
    half = provider(midpoint_file, 'midpoint_rhs')
    accepted, audit = update(initial, midpoint, half, .001, .001, backend, shadow_audit=True)
    write(args.output/'full_factor.json', audit)
    np.savez(args.output/'endpoint.npz', **accepted)
    native = load(args.native_checkpoint)
    exact = {key: bool(np.array_equal(accepted[key], native[key])) for key in native}
    differences = {key: float(np.max(np.abs(np.asarray(accepted[key], dtype=np.longdouble)
                    - np.asarray(native[key], dtype=np.longdouble)))) for key in native}
    write(args.output/'comparison.json', dict(status='passed_samehost_original_global_midpoint'
        if all(exact.values()) else 'failed_samehost_original_global_midpoint',
        host=socket.gethostname(), exact=exact, maximum_absolute_difference=differences,
        physical_providers=2, accepted_time_steps=1,
        reference='Unchanged original global-stage API and Python midpoint update on same host/GPU arithmetic library',
        no_tolerance_change=True, long_time_stability_verified=False))
    assert all(exact.values()), (exact, differences)
    # Only successful disposable work belonging to this invocation is removed.
    shutil.rmtree(memory)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=('main', 'stage'), default='main')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--state', type=Path)
    p.add_argument('--report', type=Path)
    p.add_argument('--native-checkpoint', type=Path)
    args = p.parse_args()
    if args.mode == 'stage':
        stage(args)
    else:
        main(args)
