"""Repair missing audit import; preserve two completed numerical qualifiers."""
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/vertex_halo_pair_t2654918_ada_20261029')
SOURCE = ROOT.parent/'vertex_halo_pair_t2654918_a100_20261029'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if (ROOT/'manual_exit_status.txt').read_text().strip() != '1':
        raise ValueError('Expected terminated initial adapter execution')
    if "No module named 'audit_kepler_polar_long_gate_20261027'" not in (ROOT/'manual_runner.stderr').read_text():
        raise ValueError('Different failure; refuse automatic repair')
    for name in ['qualifier_old', 'qualifier_new']:
        if 'EXIT:0' not in (ROOT/name/'log').read_text():
            raise ValueError('Completed native qualifier missing')
        gate = json.loads((ROOT/name/'audit_endpoint_gate.json').read_text())
        if not all(gate['checks'].values()):
            raise ValueError('Native physics qualifier did not pass')
    remaining = ['control40', 'vertex40', 'vertex20', 'vertex20_reuse', 'vertex20_reuse_local_sync']
    if any((ROOT/name/'log').exists() for name in remaining) or (ROOT/'qualifier.json').exists():
        raise ValueError('Unexpected partial subsequent work; refusing replay')
    archive = ROOT/'packaging_recovery_20261029'
    archive.mkdir()
    adapter = 'audit_kepler_vertex_halo_manual_20261029.py'
    for name in [adapter, 'audit_inputs.sha256', 'staging_receipt.json',
                 'manual_runner.stdout', 'manual_runner.stderr', 'manual_exit_status.txt']:
        shutil.copy2(ROOT/name, archive/name)
    helper = 'audit_kepler_polar_long_gate_20261027.py'
    if (ROOT/helper).exists():
        raise ValueError('Missing helper unexpectedly exists')
    shutil.copy2(SOURCE/helper, ROOT/helper)
    shutil.copy2(DOCS/adapter, ROOT/adapter)
    subprocess.run([sys.executable, '-c',
        'import audit_kepler_vertex_halo_manual_20261029; import audit_kepler_vertex_halo_pair_20261029'],
        cwd=ROOT, check=True)
    script = (ROOT/'run_cases.sh').read_text()
    old = 'for name in qualifier_old qualifier_new '+ ' '.join(remaining) + '; do'
    if script.count(old) != 1:
        raise ValueError('Unexpected original case sequence')
    script = script.replace(old, 'for name in '+' '.join(remaining)+'; do')
    marker = 'pair_start_seconds=$SECONDS\n'
    script = script.replace(marker, marker+f'/home/kjhan/miniconda3/bin/python3 {ROOT/adapter} --qualifier > {ROOT}/qualifier_recovery.log\n', 1)
    (ROOT/'run_cases_remaining.sh').write_text(script)
    wrapper = (ROOT/'run_manual.sh').read_text()
    marker = "export KEPLER_MANUAL_SESSION_ID=$(/home/kjhan/miniconda3/bin/python3 -c 'import uuid;print(uuid.uuid4())')"
    if wrapper.count(marker) != 1:
        raise ValueError('Unexpected manual session initializer')
    replacement = "export KEPLER_MANUAL_SESSION_ID=$(/home/kjhan/miniconda3/bin/python3 -c 'import json;print(json.load(open(\"manual_execution.json\"))[\"manual_session_id\"])')"
    wrapper = wrapper.replace(marker, replacement).replace('--record-execution', '--record-continuation')
    wrapper = wrapper.replace('exec bash run_cases.sh', 'exec bash run_cases_remaining.sh')
    (ROOT/'run_manual_remaining.sh').write_text(wrapper)
    for name in ['run_cases_remaining.sh', 'run_manual_remaining.sh']:
        subprocess.run(['bash', '-n', str(ROOT/name)], check=True)
    receipt = json.loads((ROOT/'staging_receipt.json').read_text())
    files = list(receipt['source_hashes'])+[helper, 'run_cases_remaining.sh', 'run_manual_remaining.sh']
    receipt['source_hashes'] = {name: sha(ROOT/name) for name in files}
    recovery = {'observed_at': datetime.datetime.now().astimezone().isoformat(),
                'cause': 'Missing root-level transitive audit import',
                'qualifier_native_runs_reused': True, 'completed_native_cases_replayed': 0,
                'remaining_native_cases': remaining,
                'physical_gates_changed': False,
                'execution_context': 'Same actual host/GPU and benchmark UUID; new actual supervisor recorded',
                'archive': str(archive), 'source_hashes': receipt['source_hashes']}
    receipt['packaging_recovery'] = recovery
    (ROOT/'audit_inputs.sha256').write_text(''.join(f'{sha(ROOT/name)}  {name}\n' for name in files))
    (ROOT/'staging_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (ROOT/'packaging_recovery.json').write_text(json.dumps(recovery, indent=2)+'\n')
    (DOCS/'KEPLER_VERTEX_HALO_ADA_STAGE_2026-10-29.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (DOCS/'KEPLER_VERTEX_HALO_ADA_RECOVERY_2026-10-29.json').write_text(json.dumps(recovery, indent=2)+'\n')
    print(json.dumps(recovery, indent=2))


if __name__ == '__main__':
    main()
