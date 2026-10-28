"""Prepare and bind a separate continuation; never submits or edits the original chain."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

CANDIDATE_SHA = 'a22e0e8ebea0eca807ab3762d95bf037184b7784218a42fd2eb1ae7a24526bfc'
BASE_CHECKER_SHA = '618a76315a1cafa0fd767fb2006134df907b4e15c0e7f8a844117d9ffc994841'
STAGES = [('to4p0', 4.), ('to5p0', 5.), ('to6p0', 6.), ('to6p5', 6.5),
          ('to7p0', 7.), ('to8p0', 8.), ('to9p0', 9.), ('to10p0', 10.)]


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


def verify_files(root, manifest):
    for name, expected in manifest.items():
        if sha(root/name) != expected:
            raise ValueError('Candidate file changed: '+name)


def stage_chain(original, root, binary):
    original, root, binary = map(Path, [original, root, binary])
    if sha(binary) != CANDIDATE_SHA:
        raise ValueError('Unqualified candidate executable')
    if sha(original/'check_radial_and_next.before_checkpoint_validation_20261029.py') != BASE_CHECKER_SHA:
        raise ValueError('Original radial continuation gate changed')
    root.mkdir()
    shutil.copy2(binary, root/'eunha2.exe')
    for helper in ['prepare_restart.py', 'audit_kepler_polar_long_gate_20261027.py',
                   'audit_kepler_native_radial_20261026.py']:
        shutil.copy2(original/helper, root/helper)
    shutil.copy2(original/'check_radial_and_next.before_checkpoint_validation_20261029.py',
                 root/'check_radial_and_next.py')
    device_source = '''import json,socket,subprocess,sys
m=json.load(open(sys.argv[1]))
if not m.get('bound'):raise SystemExit('Candidate device binding is not activated')
lines=subprocess.check_output(['nvidia-smi','--query-gpu=uuid,name','--format=csv,noheader'],text=True,timeout=15).strip().splitlines()
if len(lines)!=1:raise SystemExit('Expected one physical qualified GPU')
gpu,model=[x.strip() for x in lines[0].split(',',1)]
if (socket.gethostname().split('.')[0],gpu,model)!=(m['host'],m['GPU_UUID'],m['GPU_model']):raise SystemExit('Unqualified production device')
'''
    (root/'verify_device.py').write_text(device_source)
    atomic(root/'device_binding.json', {'bound': False})
    files = ['eunha2.exe', 'prepare_restart.py', 'audit_kepler_polar_long_gate_20261027.py',
             'audit_kepler_native_radial_20261026.py', 'check_radial_and_next.py',
             'verify_device.py', 'device_binding.json']
    for name, target in STAGES:
        old = original/name
        stage = root/name
        stage.mkdir()
        meta = json.loads((old/'stage.json').read_text())
        if meta['target_time'] != target:
            raise ValueError('Original scientific stop changed')
        meta['parent_gate'] = (meta['parent_gate'] if name == 'to4p0' else
                               meta['parent_gate'].replace(str(original), str(root)))
        if meta['next_job']:
            meta['next_job'] = meta['next_job'].replace(str(original), str(root))
        meta.update(binary_sha256=CANDIDATE_SHA,
                    purpose='Continue the same accepted trajectory with exact validated halo40 cache/alias optimizations',
                    single_change='Same physical solver; padding reuse=1, local alias=1, vertex certificate=0, halo=40')
        atomic(stage/'stage.json', meta)
        shutil.copy2(old/'params.dat', stage/'params.dat')
        (stage/'eunha2.exe').symlink_to(root/'eunha2.exe')
        (stage/'khout.000000.dat').symlink_to((old/'khout.000000.dat').resolve())
        script = (old/'job.slurm').read_text().replace(str(original), str(root))
        if script.count('GFS_GPU_LOCAL_PERIODIC_SYNC=0') != 1:
            raise ValueError('Unexpected original alias policy')
        script = script.replace('GFS_GPU_LOCAL_PERIODIC_SYNC=0', 'GFS_GPU_LOCAL_PERIODIC_SYNC=1')
        script = script.replace('#SBATCH --nodes=1\n', '#SBATCH --nodes=1\n#SBATCH --nodelist=syn103\n', 1)
        script = script.replace('PY=/home/kjhan/miniconda3/bin/python3\n',
            f'PY=/home/kjhan/miniconda3/bin/python3\n(cd {root} && sha256sum -c chain_files.sha256)\n'
            f'"$PY" {root}/verify_device.py {root}/device_binding.json\n', 1)
        script = script.replace('TIMEFORMAT=', 'export GFS_REUSE_TESS_PADDING=1 GFS_GPU_VERTEX_HALO_CERT=0\nTIMEFORMAT=', 1)
        (stage/'job.slurm').write_text(script)
        subprocess.run(['bash', '-n', str(stage/'job.slurm')], check=True)
        (stage/'inputs.sha256').write_text(''.join(f'{sha(stage/item)}  {item}\n'
            for item in ['eunha2.exe', 'params.dat', 'khout.000000.dat']))
        files += [name+'/'+item for item in ['stage.json', 'params.dat', 'job.slurm', 'inputs.sha256']]
    manifest = {name: sha(root/name) for name in files}
    atomic(root/'template_manifest.json', {'files': manifest, 'bound': False,
        'original_first_next_job': str(original/'to4p0/job.slurm'),
        'candidate_first_next_job': str(root/'to4p0/job.slurm'),
        'stops': [t for _, t in STAGES], 'physics_or_15_plus_9_gates_relaxed': False})
    return manifest


def bind_chain(root, template, checkpoint, step, checksum, device):
    root = Path(root)
    verify_files(root, template)
    if sha(checkpoint) != checksum:
        raise ValueError('Latest accepted checkpoint changed')
    for name, _ in STAGES:
        if any((root/name/item).exists() for item in ['log', 'restart_receipt.json', 'audit_gate.json']):
            raise ValueError('Candidate continuation already has output')
    path = root/'to4p0/stage.json'
    meta = json.loads(path.read_text())
    if meta['initial_time'] != 3.:
        raise ValueError('Candidate is not the original t3 continuation')
    meta.update(checkpoint=str(checkpoint), checkpoint_step=step, checkpoint_sha256=checksum)
    atomic(path, meta)
    atomic(root/'device_binding.json', {'bound': True, **device})
    names = list(template)
    final = {name: sha(root/name) for name in names}
    (root/'chain_files.sha256').write_text(''.join(f'{value}  {name}\n' for name, value in final.items()))
    atomic(root/'bound_manifest.json', {'files': final, 'checkpoint_sha256': checksum,
        'checkpoint': str(checkpoint), 'checkpoint_step': step,
        'chain_files_sha256': sha(root/'chain_files.sha256'), 'device': device})
    verify_files(root, final)
    return sha(root/'bound_manifest.json')
