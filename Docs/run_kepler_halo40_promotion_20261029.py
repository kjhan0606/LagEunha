"""Finite post-t3 A100 validation; select a bound successor only after all gates pass."""
import argparse
import datetime
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import socket
import struct
import subprocess
import sys
import time

FLAGS = ['GFS_GPU_VERTEX_HALO_CERT', 'GFS_KEPLER_HALO_CELLS',
         'GFS_REUSE_TESS_PADDING', 'GFS_GPU_LOCAL_PERIODIC_SYNC']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic(path, value):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n')
    tmp.replace(path)


def snapshot(path):
    path = Path(path)
    with path.open('rb') as stream:
        layout = struct.unpack('<4i', stream.read(16))
        stream.seek(-52, 2)
        footer = struct.unpack('<ddfiiddff', stream.read(52))
    if layout != (-1, 147456, 576, 1) or path.stat().st_size != 16+147456*576+52:
        raise ValueError('Incomplete native checkpoint')
    if not all(math.isfinite(x) for x in footer):
        raise ValueError('Nonfinite native footer')
    return {'time': footer[0], 'dt': footer[1], 'dtold': footer[2]}


def device():
    lines = subprocess.check_output(['nvidia-smi', '--query-gpu=uuid,name', '--format=csv,noheader'],
                                    text=True, timeout=15).strip().splitlines()
    if len(lines) != 1:
        raise ValueError('Only the one-physical-GPU production node is qualified')
    gpu, model = [x.strip() for x in lines[0].split(',', 1)]
    if not gpu.startswith('GPU-') or 'A100' not in model:
        raise ValueError('Expected actual A100 GPU')
    return {'host': socket.gethostname().split('.')[0], 'GPU_UUID': gpu, 'GPU_model': model}


def record():
    row = {**device(), **{key: os.environ.get(key, '') for key in
        FLAGS+['SLURM_JOB_ID', 'SLURM_NODELIST', 'CUDA_VISIBLE_DEVICES', 'OMP_NUM_THREADS',
               'GFS_EXACT_STOP', 'GFS_DIAG_DT_MAX', 'GFS_DIAG_CHECKPOINT_STEP']}}
    with Path('launch_context.json').open('x') as stream:
        json.dump(row, stream, indent=2)
        stream.write('\n')


def bridge_view(step):
    output = Path(f'khout.{step:06d}.dat')
    snapshot(output)
    view = Path('audit_view')
    view.mkdir()
    for name in ['khout.000000.dat', output.name, 'log']:
        (view/name).symlink_to((Path.cwd()/name).resolve())


def make_case(root, name, kind, checkpoint, step, target, dumps, bridge=False):
    """Only case-local files are created; every input and wrapper is hashed."""
    root = Path(root)
    run = root/name
    run.mkdir()
    inputs = {'eunha2.exe': root/('canonical.exe' if kind == 'canonical' else 'candidate.exe'),
              'params.dat': root/'params.dat', 'khout.000000.dat': root/'initial.dat',
              f'khout.{step:06d}.dat': Path(checkpoint)}
    for helper in ['audit_kepler_polar_long_gate_20261027.py', 'audit_kepler_native_radial_20261026.py']:
        inputs[helper] = root/helper
    for item, source in inputs.items():
        (run/item).symlink_to(source.resolve())
    script = (root/'templates'/('canonical.sh' if kind == 'canonical' else 'control.sh')).read_text()
    replacements = {'GFS_KEPLER_HALO_CELLS': '40', 'GFS_GPU_VERTEX_HALO_CERT': '0',
        'GFS_REUSE_TESS_PADDING': '1' if kind == 'combo' else '0',
        'GFS_GPU_LOCAL_PERIODIC_SYNC': '1' if kind == 'combo' else '0',
        'GFS_GEOM_FORCE_DUMP_MAXCALL': str(dumps), 'HYDRO_TSTOP': repr(target),
        'EUNHA_DUMP_DT': '1000000000' if bridge else '.05', 'GFS_EXACT_STOP': '0' if bridge else '1'}
    for key, value in replacements.items():
        script, count = re.subn(r'\b'+key+r'=\S+', key+'='+value, script)
        if count != 1:
            raise ValueError('Unexpected case template key: '+key)
    script, count = re.subn(r'^"\$PY" .*record_kepler_vertex_case_context_20261029.py$',
        ('unset GFS_DIAG_DT_MAX\n'+(f'export GFS_DIAG_CHECKPOINT_STEP={step+1}\n' if bridge else
                                   'unset GFS_DIAG_CHECKPOINT_STEP\n'))+
        f'"$PY" {root/Path(__file__).name} --record', script, flags=re.M)
    if count != 1:
        raise ValueError('Missing template context recorder')
    script, count = re.subn(r'params.dat 23008 > log', f'params.dat {step} > log', script)
    if count != 1:
        raise ValueError('Missing restart step in template')
    time_value = '"$actual_time"' if bridge else repr(target)
    script, count = re.subn(r'--time \S+ --output audit_endpoint_gate.json',
                           '--time '+time_value+' --output audit_endpoint_gate.json', script)
    if count != 1:
        raise ValueError('Missing endpoint audit target')
    if bridge:
        marker = '"$PY" ./audit_kepler_polar_long_gate_20261027.py'
        if script.count(marker) != 1:
            raise ValueError('Missing native bridge audit')
        prelude = (f'"$PY" {root/Path(__file__).name} --bridge-view {step+1}\n'
                   f'actual_time=$("$PY" {root/Path(__file__).name} --time-of khout.{step+1:06d}.dat)\n')
        script = script.replace(marker, prelude+marker, 1)
        script = script.replace('"$PWD" --time', '"$PWD/audit_view" --time')
    (run/'run.sh').write_text(script)
    subprocess.run(['bash', '-n', str(run/'run.sh')], check=True)
    (run/'inputs.sha256').write_text(''.join(f'{sha(source)}  {item}\n' for item, source in inputs.items())
                                   +f'{sha(run/"run.sh")}  run.sh\n')
    return run


def case_context(root, name, parent, gpu, expected, bridge=False):
    row = json.loads((root/name/'launch_context.json').read_text())
    checks = [row.get('SLURM_JOB_ID') == parent, row.get('OMP_NUM_THREADS') == '2',
              all(row.get(key) == value for key, value in gpu.items()),
              all(row.get(key) == str(value) for key, value in zip(FLAGS, expected))]
    if bridge:
        checks += [row.get('GFS_EXACT_STOP') == '0', not row.get('GFS_DIAG_DT_MAX')]
    return all(checks)


def branches_executed(root, name):
    text = (root/name/'log').read_text()
    return (bool(re.search(r'\[TESS_PADDING_REUSE\] rank=0 .*reuse=1 refresh_validation=1', text))
            and '[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh' in text
            and '[HALO_SYNC]' not in text and 'alias invariants failed' not in text)


def duplicate_result(report, config_hash, current_next, original_next, candidate_next):
    if (not report.get('finished_at') or not report.get('cleanup_verified')
            or not report.get('safe_to_advance')
            or report.get('config_sha256') != config_hash):
        return 32
    if report.get('status') == 'ORIGINAL_NEXT_RETAINED' and current_next == original_next:
        return 0
    if report.get('status') == 'PROMOTED' and report.get('promoted') and current_next == candidate_next:
        return 0
    return 32


def promotion_transaction(stage, before, next_job, eligibility, chain_manifest):
    """One CAS under a lock; never submits, reverts operator choices, or edits running code."""
    required = ['halo40', 'performance', 'bridge', 'checkpoint_sha256']
    if not all(re.fullmatch(r'[0-9a-f]{64}', eligibility.get(k, '')) for k in required):
        raise ValueError('Missing immutable eligibility receipts')
    if not re.fullmatch(r'[0-9a-f]{64}', chain_manifest):
        raise ValueError('Missing bound candidate-chain manifest')
    with (stage/'halo40_next_job.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        current = (stage/'stage.json').read_bytes()
        if current != before:
            raise ValueError('Original next-job metadata changed during validation')
        receipt = {'status': 'PREPARED', 'eligibility': eligibility,
                   'bound_chain_manifest_sha256': chain_manifest,
                   'old_stage_sha256': hashlib.sha256(before).hexdigest(), 'next_job': str(next_job)}
        atomic(stage/'halo40_promotion_transaction.json', receipt)
        meta = json.loads(current)
        meta['next_job'] = str(next_job)
        atomic(stage/'stage.json', meta)
        if json.loads((stage/'stage.json').read_text())['next_job'] != str(next_job):
            raise RuntimeError('Promotion readback failed')
        receipt.update(status='COMMITTED', new_stage_sha256=sha(stage/'stage.json'))
        atomic(stage/'halo40_promotion_transaction.json', receipt)
        return receipt


def run(config_path):
    import run_kepler_checkpoint_validation_20261029 as guard
    import audit_kepler_padding_reuse_pair_20261028 as common
    from audit_kepler_geometry_force_dump_20261029 import compare_dumps
    from audit_kepler_halo40_net_wall_20261029 import performance_summary
    from audit_kepler_restart_bytes_20261029 import compare as restart_compare
    from kepler_halo40_candidate_chain_20261029 import bind_chain
    config = json.loads(config_path.read_text())
    stage, root = Path(config['stage']), Path(config['root'])
    parent = str(config['parent_job'])
    if Path.cwd().resolve() != stage or os.environ.get('SLURM_JOB_ID') != parent:
        return 0
    report_path = stage/'halo40_post_checkpoint.json'
    config_hash = sha(config_path)
    try:
        with (stage/'halo40_post_checkpoint.lock').open('x') as stream:
            stream.write(str(os.getpid())+'\n')
    except FileExistsError:
        if report_path.exists():
            old = json.loads(report_path.read_text())
            return duplicate_result(old, config_hash, json.loads((stage/'stage.json').read_text())['next_job'],
                config['expected_next_job'], str(Path(config['candidate_root'])/'to4p0/job.slurm'))
        return 32
    if report_path.exists():
        return 32
    report = {'parent_job': parent, 'config_sha256': config_hash,
              'started_at': datetime.datetime.now().astimezone().isoformat(),
              'status': 'PREFLIGHT', 'promoted': False}
    before = (stage/'stage.json').read_bytes()
    process = None
    clean = True
    must_stop = False
    commit_started = False
    checkpoint = None
    checkpoint_hash = None
    try:
        gate = json.loads((stage/'audit_gate.json').read_text())
        radial = json.loads((stage/'continuation_gate.json').read_text())
        guard.gate_policy(gate, radial, os.environ, parent)
        if json.loads(before)['next_job'] != config['expected_next_job']:
            raise ValueError('Original next job was changed')
        prior = json.loads((stage/'post_checkpoint_validation.json').read_text())
        if prior.get('status') == 'RUNNING_AUXILIARY' or prior.get('cleanup_error'):
            clean = False
            raise RuntimeError('Original auxiliary did not finish cleanly')
        if prior.get('auxiliary_process_group') and not prior.get('auxiliary_steps_reaped'):
            clean = False
            raise RuntimeError('Original auxiliary numeric steps are unverified')
        if guard.steps(parent):
            clean = False
            raise RuntimeError('Unexpected original numeric MPI steps remain')
        info = guard.fields(guard.command('scontrol', 'show', 'job', parent, '-o'))
        owner = re.search(r'\((\d+)\)', info.get('UserId', ''))
        if not owner or int(owner[1]) != os.getuid():
            must_stop = True
            raise RuntimeError('Parent allocation owner mismatch')
        available = guard.remaining_seconds(os.environ, info, time.time())
        if info.get('JobState') != 'RUNNING':
            clean = False
            raise RuntimeError('Parent allocation is no longer running')
        if available < 3600:
            raise ValueError('Need at least one allocation hour remaining')
        if os.environ.get('GFS_DIAG_DT_MAX'):
            raise ValueError('Unexpected production diagnostic timestep cap')
        gpu = device()
        if gpu['host'] != 'syn103':
            raise ValueError('Unqualified production node')
        apps = subprocess.check_output(['nvidia-smi', '--query-compute-apps=pid', '--format=csv,noheader'],
                                       text=True, timeout=15).strip()
        if apps:
            clean = False
            raise RuntimeError('Unexpected GPU compute process before auxiliary start')
        for name, expected in config['frozen_hashes'].items():
            if sha(root/name) != expected:
                raise ValueError('Frozen validation source changed: '+name)
        if any((root/name/'log').exists() for name in config['static_cases']):
            raise ValueError('Partial validation output exists; do not replay')
        checkpoint = Path(gate['endpoint'])
        if snapshot(checkpoint)['time'] != 3.:
            raise ValueError('Not the original exact accepted t3')
        step = int(re.fullmatch(r'khout\.(\d+)\.dat', checkpoint.name)[1])
        checkpoint_hash = sha(checkpoint)
        deadline = time.time()+min(3000, available-600)
        common.ROOT = root
        common.RANKS.update({name: 1 for name in config['static_cases']+['reference40',
                            'bridge_canonical', 'bridge_off', 'bridge_combo']})

        def run_case(name, budget=900):
            nonlocal process
            if deadline-time.time() < budget+240:
                raise TimeoutError('Insufficient case and cleanup budget; interval is not shortened')
            case = root/name
            report.update(status='RUNNING_CASE', current_case=name)
            atomic(report_path, report)
            with (case/'runner.out').open('xb') as out, (case/'runner.err').open('xb') as err:
                process = subprocess.Popen(['bash', 'run.sh'], cwd=case, stdout=out, stderr=err,
                                           start_new_session=True)
                try:
                    code = process.wait(timeout=min(budget, deadline-time.time()-240))
                except subprocess.TimeoutExpired:
                    guard.stop_group(process)
                    raise TimeoutError('Native case timeout: '+name)
            if code:
                raise RuntimeError('Native case or physics auditor failed: '+name)

        cached = Path(config['cached_control'])
        can_reuse = False
        try:
            context = json.loads((cached/'launch_context.json').read_text())
            can_reuse = (context['SLURM_JOB_ID'] == parent and context['SLURM_NODELIST'] == 'syn103'
                         and context.get('CUDA_VISIBLE_DEVICES', '') == os.environ.get('CUDA_VISIBLE_DEVICES', '')
                         and all(json.loads((cached/'audit_endpoint_gate.json').read_text())['checks'].values())
                         and all(context.get(k) == str(v) for k, v in zip(FLAGS, (0,40,0,0)))
                         and sha(cached/'eunha2.exe') == config['candidate_binary_sha256']
                         and all((cached/f'geom_force_rank000_call{i:03d}.bin').is_file() for i in range(1, 5)))
            if can_reuse:
                subprocess.run(['sha256sum', '-c', 'inputs.sha256'], cwd=cached,
                               capture_output=True, check=True, timeout=60)
        except (OSError, ValueError, KeyError):
            pass
        if not can_reuse:
            run_case('fresh_reference40')
            cached = root/'fresh_reference40'
        (root/'reference40').symlink_to(cached)
        report['cached_reference_reused'] = can_reuse
        run_case('combined40')
        combo_state = common.compare('reference40', 'combined40', exact=True)
        combo_geometry = [compare_dumps(root/'reference40'/f'geom_force_rank000_call{i:03d}.bin',
                          root/'combined40'/f'geom_force_rank000_call{i:03d}.bin') for i in range(1, 5)]
        combo_pass = combo_state['passed'] and all(x['passed'] for x in combo_geometry)
        combo_pass &= case_context(root, 'combined40', parent, gpu, (0, 40, 1, 1))
        combo_pass &= branches_executed(root, 'combined40')
        atomic(root/'halo40_comparison.json', {'passed': bool(combo_pass), 'state': combo_state, 'geometry': combo_geometry})
        if not combo_pass:
            raise ValueError('Exact A100 halo40 combination gate failed')
        for name in ['old_a', 'combined_b', 'combined_b2', 'old_a2']:
            run_case(name)
        pairs = [common.compare('old_a', 'combined_b', exact=True),
                 common.compare('combined_b', 'combined_b2', exact=True),
                 common.compare('old_a', 'old_a2', exact=True)]
        walls = {}
        perf_context = True
        for name in ['old_a', 'combined_b', 'combined_b2', 'old_a2']:
            _, _, metrics = common.load(name)
            walls[name] = metrics['wall_seconds']
            perf_context &= metrics['steps'] == 59 and abs(metrics['time']-2.66) < 1e-12
            expected = (0, 40, 1, 1) if name.startswith('combined') else (0, 40, 0, 0)
            perf_context &= case_context(root, name, parent, gpu, expected)
            perf_context &= not any((root/name).glob('geom_force_rank*.bin'))
            if name.startswith('combined'):
                perf_context &= branches_executed(root, name)
        perf_endpoints = {name: Path(json.loads((root/name/'audit_endpoint_gate.json').read_text())['endpoint'])
                          for name in ['old_a', 'combined_b', 'combined_b2', 'old_a2']}
        perf_bytes = [restart_compare(perf_endpoints['old_a'], perf_endpoints[name])
                      for name in ['combined_b', 'combined_b2', 'old_a2']]
        performance = performance_summary(walls)
        performance['scope'] = 'This actual A100 allocation and interval; not a long-time stability claim'
        perf_pass = (all(p['passed'] for p in pairs+perf_bytes) and perf_context
                     and performance['clear_improvement'])
        atomic(root/'net_wall_comparison.json', {'passed': bool(perf_pass), 'performance': performance,
                                               'pairs': pairs, 'full_restart_bytes': perf_bytes})
        if not perf_pass:
            raise ValueError('Canonical A100 wall improvement is not established')
        bridge_target = math.nextafter(3., math.inf)
        for name, kind in [('bridge_canonical', 'canonical'), ('bridge_off', 'off'), ('bridge_combo', 'combo')]:
            make_case(root, name, kind, checkpoint, step, bridge_target, 0 if kind == 'canonical' else 4, bridge=True)
            run_case(name, 180)
        bridge_pairs = [common.compare('bridge_canonical', 'bridge_off', exact=True),
                        common.compare('bridge_off', 'bridge_combo', exact=True)]
        bridge_geometry = [compare_dumps(root/'bridge_off'/f'geom_force_rank000_call{i:03d}.bin',
                            root/'bridge_combo'/f'geom_force_rank000_call{i:03d}.bin') for i in range(1, 5)]
        endpoints = {name: Path(json.loads((root/name/'audit_endpoint_gate.json').read_text())['endpoint'])
                     for name in ['bridge_canonical', 'bridge_off', 'bridge_combo']}
        bridge_bytes = [restart_compare(endpoints['bridge_canonical'], endpoints[name])
                        for name in ['bridge_off', 'bridge_combo']]
        bridge_pass = all(p['passed'] for p in bridge_pairs+bridge_geometry+bridge_bytes)
        for name, expected in [('bridge_canonical', (0,40,0,0)), ('bridge_off', (0,40,0,0)), ('bridge_combo', (0,40,1,1))]:
            _, trace, metrics = common.load(name)
            end = snapshot(endpoints[name])
            bridge_pass &= (metrics['steps'] == 1 and int(trace['steps'][0][1]) == step+1
                            and 3. < end['time'] < 4. and end['dt'] > bridge_target-3.)
            bridge_pass &= case_context(root, name, parent, gpu, expected, bridge=True)
        bridge_pass &= branches_executed(root, 'bridge_combo')
        atomic(root/'bridge_comparison.json', {'passed': bool(bridge_pass), 'pairs': bridge_pairs,
               'geometry': bridge_geometry, 'full_restart_bytes': bridge_bytes,
               'production_checkpoint': str(checkpoint), 'production_checkpoint_sha256': checkpoint_hash,
               'method': 'EXACT_STOP=0, no diagnostic cap, nextafter loop guard and explicit accepted-step dump; natural production CFL step'})
        if not bridge_pass:
            raise ValueError('Actual next production-step bridge failed')
        clean, actions = guard.cleanup_steps(parent)
        report['cleanup_actions'] = actions
        if not clean:
            raise RuntimeError('Unverified residual auxiliary ranks')
        if time.time() > deadline:
            raise TimeoutError('Auxiliary validation deadline expired')
        if sha(checkpoint) != checkpoint_hash:
            must_stop = True
            raise ValueError('Accepted checkpoint changed before promotion')
        for name, expected in config['frozen_hashes'].items():
            if sha(root/name) != expected:
                raise ValueError('Validation source changed during execution: '+name)
        for name in ['halo40_comparison.json', 'net_wall_comparison.json', 'bridge_comparison.json']:
            if json.loads((root/name).read_text()).get('passed') is not True:
                raise ValueError('Eligibility receipt is missing, stale or failed: '+name)
        eligibility = {'halo40': sha(root/'halo40_comparison.json'),
                       'performance': sha(root/'net_wall_comparison.json'),
                       'bridge': sha(root/'bridge_comparison.json'),
                       'checkpoint_sha256': checkpoint_hash, 'device': gpu}
        chain_hash = bind_chain(Path(config['candidate_root']), config['candidate_template'],
                                checkpoint, step, checkpoint_hash, gpu)
        if time.time() > deadline:
            raise TimeoutError('Promotion deadline expired during chain binding')
        commit_started = True
        transition = promotion_transaction(stage, before, Path(config['candidate_root'])/'to4p0/job.slurm',
                                           eligibility, chain_hash)
        report.update(status='PROMOTED', promoted=True, transition=transition)
    except BaseException as error:
        report['error'] = repr(error)
        report['status'] = 'OPERATIONAL_STOP' if commit_started or must_stop else 'ORIGINAL_NEXT_RETAINED'
        if process is not None:
            try:
                if process.poll() is None:
                    guard.stop_group(process)
                clean, actions = guard.cleanup_steps(parent)
                report['cleanup_actions'] = actions
            except BaseException as cleanup_error:
                clean = False
                report['cleanup_error'] = repr(cleanup_error)
        if isinstance(error, (KeyboardInterrupt, InterruptedError, SystemExit)):
            must_stop = True
    finally:
        if checkpoint is not None and checkpoint_hash is not None:
            try:
                intact = sha(checkpoint) == checkpoint_hash
            except OSError:
                intact = False
            report['checkpoint_integrity_verified'] = intact
            if not intact:
                must_stop = True
                report['status'] = 'OPERATIONAL_STOP'
        report['cleanup_verified'] = clean
        report['safe_to_advance'] = clean and not must_stop and (not commit_started or report['promoted'])
        if not report['safe_to_advance']:
            report['status'] = 'OPERATIONAL_STOP'
        report['finished_at'] = datetime.datetime.now().astimezone().isoformat()
        atomic(report_path, report)
    return 0 if report['safe_to_advance'] else 32


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--run', type=Path)
    modes.add_argument('--record', action='store_true')
    modes.add_argument('--time-of', type=Path)
    modes.add_argument('--bridge-view', type=int)
    args = parser.parse_args()
    if args.record:
        record()
    elif args.time_of:
        print(snapshot(args.time_of)['time'])
    elif args.bridge_view is not None:
        bridge_view(args.bridge_view)
    else:
        def interrupted(signum, frame):
            raise InterruptedError('Signal '+str(signum))
        signal.signal(signal.SIGTERM, interrupted)
        signal.signal(signal.SIGINT, interrupted)
        raise SystemExit(run(args.run))
