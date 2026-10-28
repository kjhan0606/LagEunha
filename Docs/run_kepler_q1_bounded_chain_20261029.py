"""Finite checkpoint chain within one Slurm allocation; never a monitor loop."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n');temporary.replace(path)


def run(root):
    root = root.resolve()
    cfg = json.loads((root/'chain.json').read_text())
    if os.environ.get('SLURM_JOB_ID') is None or Path.cwd().resolve() != root:
        raise ValueError('Run only inside the intended Slurm working directory')
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0')) != 8:
        raise ValueError('The qualified request uses eight OpenMP threads')
    if not json.loads(Path(cfg['preflight']).read_text())['passed']:
        raise ValueError('Native initial-state preflight did not pass')
    if not json.loads(Path(cfg['ID_regression']).read_text())['passed']:
        raise ValueError('Native gapped-ID collision regression did not pass')
    for name, expected in cfg['hashes'].items():
        if digest(Path(name)) != expected:
            raise ValueError('Frozen input/source changed: '+name)
    state = {'job': os.environ['SLURM_JOB_ID'], 'status': 'running', 'completed': [],
             'initial': cfg['initial'], 'policy': 'Diagnostic checkpoints .825/2.5; unchanged 15+9 physics gates mandatory from exact t=3.'}
    report = root/'chain_state.json'
    # A second dispatch must not turn an existing run into a t=0 replay.
    with (root/'execution_claim.json').open('x') as stream:
        json.dump({'job': state['job'], 'started_unix': time.time()}, stream)
    atomic(report,state)
    initial = Path(cfg['initial']);previous = initial;step = 1
    try:
        for item in cfg['stages']:
            stage = root/item['directory']
            if (stage/'log').exists() or (stage/'checkpoint_gate.json').exists():
                raise ValueError('Existing stage output: never overwrite or replay')
            checkpoint = stage/f'khout.{step:06d}.dat'
            if checkpoint.exists() or checkpoint.is_symlink():
                raise ValueError('Existing input link; explicit resume audit required')
            checkpoint.symlink_to(previous)
            atomic(stage/'restart_receipt.json', {'source': str(previous),'sha256':digest(previous),'step':step,'target_time':item['time']})
            state.update(active_stage=item['directory'],restart=str(previous));atomic(report,state)
            started = time.monotonic()
            native = subprocess.run(['bash', str(stage/'run.sh'), str(step)],cwd=stage)
            if native.returncode:
                raise RuntimeError('Native solver failed; accepted checkpoints preserved')
            audit = subprocess.run([sys.executable, str(root/'audit_kepler_q1_bounded_gate_20261029.py'),
                                    str(initial),str(stage),'--time',str(item['time']),
                                    '--output',str(stage/'checkpoint_gate.json')],cwd=stage,
                                   stdout=(stage/'checkpoint_gate.log').open('w'),stderr=subprocess.STDOUT)
            gate = json.loads((stage/'checkpoint_gate.json').read_text())
            physical = gate['checks']
            # Early output is a comparison checkpoint, not a replacement for
            # the decisive exact t=3 gate. Invalid states never continue.
            hard = ['finite','mass_exact','W0_exact','positive_U_and_volume',
                    'rho_equals_mass_over_volume','no_floor_or_dual_energy_injection',
                    'endpoint_area_tiles_polygon','disk_population_preserved',
                    'all_endpoints_inside_polygon_to_roundoff','no_unresolved_event_or_rollback_failure']
            if any(not physical[k] for k in hard) or (item['time'] >= 3. and audit.returncode != 0):
                raise RuntimeError('Checkpoint physics gate failed; do not advance or loosen thresholds')
            previous = Path(gate['endpoint']);step = int(previous.name.split('.')[1])
            state['completed'].append({'time':gate['time'],'endpoint':str(previous),'sha256':digest(previous),
                                       'all_checks_passed':all(physical.values()),'elapsed_seconds':time.monotonic()-started})
            atomic(report,state)
        state['status']='finite_chain_complete'
    except BaseException as exc:
        state.update(status='stopped',error=str(exc));atomic(report,state);raise
    atomic(report,state)


if __name__ == '__main__':
    run(Path(sys.argv[1]))
