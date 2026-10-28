#!/usr/bin/env python3
"""Prepare one immutable Slurm segment from a provenance-verified checkpoint.

This creates configuration only; it does not submit a job or advance state.
The generated batch job runs the native MPI/CUDA controller, audits the final
accepted state, then seals it for the next segment.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import sys


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require_token_path(path):
    value = str(Path(path).resolve())
    if any(c.isspace() for c in value) or len(value) >= 500:
        raise ValueError(f'Native parameter path contains whitespace or is too long: {value}')
    return value


def write_new(path, content, mode=0o600):
    path = Path(path)
    with path.open('x') as f:
        f.write(content)
        f.flush()
        os.fsync(f.fileno())
    path.chmod(mode)


def verify_provider(provider_dir, executable, operator, operator_sha):
    provenance_path = provider_dir / 'provenance.json'
    provenance = json.loads(provenance_path.read_text())
    sources = provenance['sources']
    for name, digest in sources.items():
        if sha(provider_dir / name) != digest:
            raise ValueError(f'Provider source changed: {name}')
    required = {
        'full_moment_provider_service.py',
        'full_moment_restart_provenance.py',
        'run_full_moment_with_identity.py',
        'liblag_straight_faces.so',
        'liblag_transport_root.so',
    }
    if not required.issubset(sources):
        raise ValueError('Provider snapshot is missing restart/runtime identity modules')
    binary_sha = sha(executable)
    if binary_sha != provenance['binary_sha256']:
        raise ValueError('Executable does not match provider snapshot provenance')
    actual_operator_sha = sha(operator)
    if actual_operator_sha != operator_sha:
        raise ValueError('CUDA operator does not match requested SHA')
    return provenance, binary_sha, sha(provenance_path)


def render_params(*, python, provider, checkpoint, initial, output, physics,
                  physics_sha, workers, threads, h, steps, restart_manifest):
    rows = [
        'define Simulation Model = FullMomentKepler',
        'moment_native_evolution = 1',
        f'moment_python = {python}',
        f'moment_provider_script = {provider}',
        f'moment_state = {checkpoint}',
        f'moment_initial = {initial}',
        f'moment_output = {output}',
        f'moment_physics_file = {physics}',
        f'moment_physics_sha = {physics_sha}',
        'moment_integrator = explicit_midpoint',
        f'moment_workers = {workers}',
        f'moment_threads = {threads}',
        f'moment_h = {h:.17g}',
        f'moment_steps = {steps}',
        'moment_validation = 0',
        f'moment_restart_manifest = {restart_manifest}',
    ]
    return '\n'.join(rows) + '\n'


def render_submit(*, run_dir, runtime, python, provider_dir, executable,
                  initial, physics, operator, operator_sha, binary_sha,
                  target_step, target_time, ranks, workers, threads, job_name):
    cpus = 8
    memory = 128 * ranks
    launcher = provider_dir / 'run_full_moment_with_identity.py'
    provider = provider_dir / 'full_moment_provider_service.py'
    restart_helper = provider_dir / 'full_moment_restart_provenance.py'
    runtime_log = runtime / 'native.log'
    identity = run_dir / 'runtime_identity.json'
    audit = runtime / f'completed_{target_step:04d}.json'
    checkpoint = runtime / f'accepted_{target_step:04d}'
    seal_dir = run_dir / 'restart_seal'

    q = shlex.quote
    return f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks={ranks}
#SBATCH --cpus-per-task={cpus}
#SBATCH --gres=gpu:H200:{ranks}
#SBATCH --mem={memory}G
#SBATCH --time=1-23:59:00
#SBATCH --chdir={run_dir}
#SBATCH --output={run_dir}/slurm-%j.out
#SBATCH --error={run_dir}/slurm-%j.out
set -euo pipefail
module purge
module load intel/tbb intel/compiler-rt intel/umf intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2
export LD_LIBRARY_PATH=/home/kjhan/local/lib:/opt/ohpc/pub/compiler/gcc/13.2.0/lib64:${{LD_LIBRARY_PATH:-}}
export OMP_NUM_THREADS={threads} OMP_PROC_BIND=false OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 I_MPI_PIN_DOMAIN=omp
export EUNHA_MOMENT_BINARY_SHA={binary_sha}
export EUNHA_MOMENT_OPERATOR_PATH={operator}
export EUNHA_MOMENT_OPERATOR_SHA={operator_sha}
export LAGFORCE_NATIVE_FACE_STAGE=1
export LAGFORCE_NATIVE_STRAIGHT_FACE_LIB={q(str(provider_dir / 'liblag_straight_faces.so'))}
export LAGFORCE_NATIVE_TRANSPORT_ROOT_LIB={q(str(provider_dir / 'liblag_transport_root.so'))}
{q(str(python))} {q(str(launcher))} \\
  --exe {q(str(executable))} --provider {q(str(provider))} \\
  --params {q(str(run_dir / 'params.dat'))} --ranks {ranks} \\
  --log {q(str(runtime_log))} --identity {q(str(identity))} \\
  > {q(str(run_dir / 'launch.log'))} 2>&1
{q(str(python))} - {q(str(runtime))} {target_step} {target_time:.17g} {ranks} <<'PY'
import json, math, sys
from pathlib import Path
runtime=Path(sys.argv[1]); step=int(sys.argv[2]); target=float(sys.argv[3]); ranks=int(sys.argv[4])
report=json.loads((runtime/f'completed_{{step:04d}}.json').read_text())
assert report['status']=='completed_native_evolution_batch'
assert report['accepted_commit'] and math.isclose(report['time'],target,rel_tol=0.0,abs_tol=32*sys.float_info.epsilon*max(1.0,abs(target))) and report['step']==step
assert report['ranks']==ranks and report['fixed_mass'] and report['Re_infinity']
assert not report['long_time_stability_verified']
assert (runtime/f'accepted_{{step:04d}}').is_dir()
PY
{q(str(python))} {q(str(restart_helper))} \\
  --checkpoint {q(str(checkpoint))} --initial {q(str(initial))} \\
  --physics {q(str(physics))} --output {q(str(seal_dir))}
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-checkpoint', type=Path, required=True)
    p.add_argument('--source-manifest', type=Path, required=True)
    p.add_argument('--initial', type=Path, required=True)
    p.add_argument('--physics', type=Path, required=True)
    p.add_argument('--provider-dir', type=Path, required=True)
    p.add_argument('--executable', type=Path, required=True)
    p.add_argument('--operator', type=Path, required=True)
    p.add_argument('--operator-sha', required=True)
    p.add_argument('--target-time', type=float, required=True)
    p.add_argument('--run-dir', type=Path, required=True)
    p.add_argument('--ranks', type=int, choices=(1, 2), default=2)
    p.add_argument('--workers', type=int, default=7)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--python', type=Path, default=Path('/home/kjhan/miniconda3/bin/python3'))
    args = p.parse_args()

    checkpoint = args.source_checkpoint.resolve()
    source_manifest = args.source_manifest.resolve()
    initial = args.initial.resolve()
    physics = args.physics.resolve()
    provider_dir = args.provider_dir.resolve()
    executable = args.executable.resolve()
    operator = args.operator.resolve()
    python = args.python.resolve()
    run_dir = args.run_dir.resolve()
    if not run_dir.parent.is_dir() or run_dir.exists():
        raise FileExistsError(f'Run directory parent must exist and destination must be new: {run_dir}')
    if not (checkpoint.is_dir() and source_manifest.is_file() and initial.is_file()
            and physics.is_file() and provider_dir.is_dir() and executable.is_file()
            and operator.is_file() and python.is_file()):
        raise FileNotFoundError('A requested source, executable, provider, operator or Python file is missing')
    if args.workers < 1 or args.threads < 1:
        raise ValueError('workers and threads must be positive')

    provider_provenance, binary_sha, provider_provenance_sha = verify_provider(
        provider_dir, executable, operator, args.operator_sha)
    os.environ['EUNHA_MOMENT_OPERATOR_PATH'] = str(operator)
    os.environ['EUNHA_MOMENT_OPERATOR_SHA'] = args.operator_sha
    sys.path.insert(0, str(provider_dir))
    import full_moment_restart_provenance as restart_provenance
    if not hasattr(restart_provenance, 'verify_source_runtime_identity'):
        raise ValueError('Provider snapshot restart helper lacks source-runtime identity verification')
    verify = restart_provenance.verify
    source_record = json.loads(source_manifest.read_text())
    required_identity = ('source_runtime_identity_path', 'source_runtime_identity_sha256',
        'source_params_path', 'source_params_sha256', 'source_binary_sha256',
        'source_operator_sha256', 'source_provider_snapshot_sources_sha256', 'ranks')
    if any(name not in source_record for name in required_identity):
        raise ValueError('Source manifest lacks a completed runtime identity; reseal it before continuation')
    h = float(source_record['h'])
    if not math.isfinite(args.target_time) or args.target_time <= float(source_record['time']):
        raise ValueError('Target time must be finite and later than the source checkpoint')
    target_step = round(args.target_time / h)
    time_tol = 32 * sys.float_info.epsilon * max(1.0, abs(args.target_time))
    if target_step < 1 or target_step > 1_000_000 or abs(target_step * h - args.target_time) > time_tol:
        raise ValueError('Target time must lie exactly on the authenticated timestep grid')
    verify(checkpoint, source_manifest, initial, physics, h)
    if target_step <= int(source_record['step']):
        raise ValueError('Target step must be later than the source checkpoint step')

    paths = {
        'python': require_token_path(python),
        'provider': require_token_path(provider_dir / 'full_moment_provider_service.py'),
        'checkpoint': require_token_path(checkpoint),
        'initial': require_token_path(initial),
        'output': require_token_path(run_dir / 'runtime'),
        'physics': require_token_path(physics),
        'restart_manifest': require_token_path(source_manifest),
    }
    params = render_params(
        python=paths['python'], provider=paths['provider'], checkpoint=paths['checkpoint'],
        initial=paths['initial'], output=paths['output'], physics=paths['physics'],
        physics_sha=source_record['physics_sha256'], workers=args.workers,
        threads=args.threads, h=h, steps=target_step,
        restart_manifest=paths['restart_manifest'])
    slug = ''.join(c if c.isalnum() or c in '-_' else '-' for c in run_dir.name)[:50]
    if not slug:
        raise ValueError('Run directory basename cannot form a Slurm job name')
    runtime = run_dir / 'runtime'
    plan = dict(
        status='prepared_native_restart_segment_not_submitted',
        source_checkpoint=str(checkpoint), source_manifest=str(source_manifest),
        source_manifest_sha256=sha(source_manifest), source_time=source_record['time'],
        source_step=source_record['step'], source_h=h,
        source_runtime_identity_path=source_record.get('source_runtime_identity_path'),
        source_runtime_identity_sha256=source_record.get('source_runtime_identity_sha256'),
        source_params_path=source_record.get('source_params_path'),
        source_params_sha256=source_record.get('source_params_sha256'),
        source_binary_sha256=source_record.get('source_binary_sha256'),
        source_operator_sha256=source_record.get('source_operator_sha256'),
        target_time=args.target_time, target_step=target_step,
        ranks=args.ranks, workers=args.workers, threads=args.threads,
        executable=str(executable), binary_sha256=binary_sha,
        provider_dir=str(provider_dir), provider_provenance_sha256=provider_provenance_sha,
        provider_source_count=len(provider_provenance['sources']),
        operator=str(operator), operator_sha256=args.operator_sha,
        native_face_stage=True, initial=str(initial), physics=str(physics),
        physics_sha256=source_record['physics_sha256'],
        accepted_checkpoint=str(runtime / f'accepted_{target_step:04d}'),
        restart_seal=str(run_dir / 'restart_seal' / 'manifest.json'),
        physical_providers=0, accepted_time_steps=0)
    submit = render_submit(
        run_dir=run_dir, runtime=runtime, python=python, provider_dir=provider_dir,
        executable=executable, initial=initial, physics=physics, operator=operator,
        operator_sha=args.operator_sha, binary_sha=binary_sha,
        target_step=target_step, target_time=args.target_time, ranks=args.ranks,
        workers=args.workers, threads=args.threads, job_name=slug)

    run_dir.mkdir(mode=0o700)
    try:
        runtime.mkdir(mode=0o700)
        write_new(run_dir / 'params.dat', params)
        write_new(run_dir / 'segment_plan.json', json.dumps(plan, indent=2) + '\n')
        write_new(run_dir / 'submit.sh', submit, 0o700)
        fd = os.open(run_dir, os.O_RDONLY | os.O_DIRECTORY)
        os.fsync(fd)
        os.close(fd)
    except Exception:
        # Keep any partial preparation visible for inspection; never overwrite
        # or silently remove a path that may contain useful evidence.
        raise
    print(json.dumps(dict(status=plan['status'], run_dir=str(run_dir),
        target_step=target_step, target_time=args.target_time,
        ranks=args.ranks, binary_sha256=binary_sha,
        operator_sha256=args.operator_sha), sort_keys=True))


if __name__ == '__main__':
    main()
