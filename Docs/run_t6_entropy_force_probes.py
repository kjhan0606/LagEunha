"""Terminal force/geometry/HG assembly from preserved matched t6 snapshots."""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import sha
b = Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--cross', action='store_true')
args=parser.parse_args()
root = b/('Codex_t6_entropy_force_cross_20261018' if args.cross else 'Codex_t6_entropy_force_probes_20261018')
binary = b/'Codex_adaptive_qr_periodic_build_20261018/eunha2.exe'
assert sha(binary) == 'e033bcc1ffbd73df73bea2b916cf204b00d3b0f9a5987ccd9cd54f1970299e7c'
assert shutil.disk_usage(b).free > 4*1024**3
root.mkdir(exist_ok=False)
template = json.loads((b/'Codex_adaptive_qr_periodic_probes_20261018/on_plus/environment.json').read_text())
manifest = {'binary': str(binary), 'binary_sha256': sha(binary), 'scope': 'Terminal assembly only; no accepted time evolution.', 'cases': {}}
for name, source, enabled in [
    ('adaptive', b/'Codex_adaptive_qr_physical_to_t10_20261018/adaptive_from_t2', 1),
    ('legacy', b/'Codex_smooth_qr_physical_to_t10_20261018/legacy_continue', 0)]:
    d = root/name
    d.mkdir()
    for filename in ('params.dat', 'khout.001783.dat'):
        shutil.copyfile(source/filename, d/filename)
    envcase = dict(template)
    if args.cross: enabled=1-enabled
    envcase.pop('GFS_QR_OPERATOR_DUMP', None)
    envcase.update(GFS_QR_SMOOTH=str(enabled), GFS_QR_ADAPTIVE=str(enabled),
                   HYDRO_TSTOP='6.001405917763673', GFS_MAX_DT='0.0025',
                   GFS_SHF_DUMP_FILE=str(d/'raw_faces.dat'))
    (d/'environment.json').write_text(json.dumps(envcase, indent=2))
    env = {k:v for k,v in os.environ.items() if not k.startswith(('GFS_', 'EUNHA_', 'SEDOV_', 'LAGMFM_', 'HYDRO_', 'OMP_'))}
    env.update(envcase)
    manifest['cases'][name] = {'source': str(source/'khout.001783.dat'), 'source_sha256': sha(d/'khout.001783.dat'), 'adaptive_operator': bool(enabled)}
    (root/'manifest.json').write_text(json.dumps(manifest, indent=2))
    with (d/'log').open('x') as stream:
        job = subprocess.Popen(['timeout','600','mpirun','-np','1',str(binary),'params.dat','1783'], cwd=d, env=env, stdout=stream, stderr=subprocess.STDOUT)
        print(name, job.pid, flush=True)
        code = job.wait()
    (d/'exit.code').write_text(str(code))
    assert code == 0, (name, code)
    assert '[ACCEPTED_STATE]' not in (d/'log').read_text()
    assert sha(d/'khout.001783.dat') == manifest['cases'][name]['source_sha256']
    print(name, 'terminal assembly complete', flush=True)
