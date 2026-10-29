"""Build a separate native binary fixing thermal64's RK4 final anchor."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def build(root):
    source = Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/project')
    baseline = source.parent/'candidate_v2/Exam/exam.c'
    if sha(baseline) != '328982249341cc46b365eb89850683d575ebe54029e453159e73a35297b3ece7':
        raise ValueError('Unexpected frozen baseline host source')
    old = 'double u_final=(direct_u?sbp[i].stress.E_inv_xx:(double)sbp[i].ie)+u_increment;'
    new = 'double u_final=(thermal64?GFS_T64_ANCHOR(sbp[i]):(direct_u?sbp[i].stress.E_inv_xx:(double)sbp[i].ie))+u_increment;'
    original = baseline.read_text()
    if original.count(old) != 1: raise ValueError('Unexpected RK4 final expression')
    dest = root/'project'
    if dest.exists(): raise FileExistsError(dest)
    dest.mkdir(parents=True)
    for p in source.iterdir():
        if p.name in ('Exam', 'eunha2.exe'): continue
        if p.name == 'eunha2.o': shutil.copy2(p, dest/p.name)
        else: (dest/p.name).symlink_to(p, target_is_directory=p.is_dir())
    exam = dest/'Exam'; exam.mkdir()
    for p in (source/'Exam').iterdir():
        if p.is_dir(): (exam/p.name).symlink_to(p, target_is_directory=True)
        elif p.is_file(): shutil.copy2(p, exam/p.name)
    (exam/'exam.c').write_text(original.replace(old, new))
    os.utime(exam/'exam.c', None)
    for name, cmd in [('exam', ['make', '-C', 'Exam', '-j2', 'USE_CUDA=1', 'GFS_VARPRESS_W0=1', 'libexam.a']),
                      ('link', ['make', '-j2', 'USE_CUDA=1', 'GFS_VARPRESS_W0=1', 'this'])]:
        with (root/('build_'+name+'.log')).open('x') as log:
            subprocess.run(cmd, cwd=dest, stdout=log, stderr=subprocess.STDOUT, check=True)
    result = {'baseline_sha256': sha(baseline), 'candidate_host_sha256': sha(exam/'exam.c'),
              'candidate_binary_sha256': sha(dest/'eunha2.exe'),
              'source_change': {'before': old, 'after': new},
              'scope': 'Experimental build only. Legacy and wall compatibility guards remain unchanged. No physical or runtime validation implied.'}
    with (root/'build_result.json').open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__': build(Path(sys.argv[1]).resolve())
