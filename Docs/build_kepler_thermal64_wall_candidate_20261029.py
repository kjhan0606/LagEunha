"""Build a separate experimental bounded thermal64 wall-transaction candidate."""
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

def replace_once(path, old, new):
    text = path.read_text()
    if text.count(old) != 1: raise ValueError('Unexpected source context: '+old[:90])
    path.write_text(text.replace(old, new))
    os.utime(path, None)

def build(root):
    source = Path('/gpfs/kjhan/LagForce/Kepler_Thermal64AnchorFix_20261029/project')
    if sha(source/'Exam/exam.c') != '09286ebd3234ae03658695f6f08a38d05958e392f83dc51f6a9c427076e2e2e4':
        raise ValueError('Unexpected anchor-fix source')
    dest = root/'project'
    if dest.exists(): raise FileExistsError(dest)
    dest.mkdir(parents=True)
    for p in source.iterdir():
        if p.name in ('Exam', 'eunha2.exe'): continue
        if p.name == 'eunha2.o': shutil.copy2(p, dest/p.name)
        else: (dest/p.name).symlink_to(p.resolve(), target_is_directory=p.is_dir())
    exam = dest/'Exam'; exam.mkdir()
    for p in (source/'Exam').iterdir():
        if p.is_dir(): (exam/p.name).symlink_to(p.resolve(), target_is_directory=True)
        elif p.is_file(): shutil.copy2(p, exam/p.name)
    c = exam/'exam.c'; wall = exam/'gfs_wall_events.inc'; thermal = exam/'gfs_thermal64_native.inc'
    replace_once(thermal, 'const char *geomface=getenv("GFS_GEOM_FACE_VEL"),*initonly=getenv("GFS_THERMAL64_INIT_ONLY");',
        'const char *geomface=getenv("GFS_GEOM_FACE_VEL"),*initonly=getenv("GFS_THERMAL64_INIT_ONLY");\n    const char *wall=getenv("GFS_KEPLER_WALL_EVENTS");')
    replace_once(thermal, '(init_gate || (tx && atoi(tx)>0 && (!txdiag || atoi(txdiag)==0)))',
        '(init_gate || (tx && atoi(tx)>0 && (!txdiag || atoi(txdiag)==0)) ||\n        (wall && strcmp(wall,"1")==0 && (!tx || atoi(tx)==0)))')
    replace_once(wall, 'a[3]+=p[i].ie+.5*m*(vx*vx+vy*vy)+m*kepler_phi(p[i].x,p[i].y);\n        a[4]+=m;a[5]+=p[i].ie;',
        'double U=gfs_t64_U(p+i);\n        a[3]+=U+.5*m*(vx*vx+vy*vy)+m*kepler_phi(p[i].x,p[i].y);\n        a[4]+=m;a[5]+=U;')
    replace_once(wall, 'gfs_vaud_on() || gfs_thermal64_requested() ||', 'gfs_vaud_on() ||')
    replace_once(wall, 'if(w->base[id].present || !isfinite(p[i].mass)',
        'if((gfs_thermal64_requested() && (!GFS_T64_TAGGED(p[i]) || !isfinite(GFS_T64_U(p[i])) || !(GFS_T64_U(p[i])>0) || !isfinite(GFS_T64_V(p[i])) || !(GFS_T64_V(p[i])>0) || !isfinite(GFS_T64_ANCHOR(p[i])))) ||\n           w->base[id].present || !isfinite(p[i].mass)')
    replace_once(wall, 'if(!(p[i].ie>0) || !(p[i].volume>0) || !isfinite(p[i].mass)',
        'if((gfs_thermal64_requested() && (!GFS_T64_TAGGED(p[i]) || !isfinite(GFS_T64_U(p[i])) || !(GFS_T64_U(p[i])>0) || !isfinite(GFS_T64_V(p[i])) || !(GFS_T64_V(p[i])>0) || !isfinite(GFS_T64_ANCHOR(p[i])) || (rhs && !isfinite(GFS_T64_RHS(p[i]))))) ||\n           !(p[i].ie>0) || !(p[i].volume>0) || !isfinite(p[i].mass)')
    replace_once(c, '    if(thermal64 && !tx.on){ fprintf(stderr,"[THERMAL64] transactional RK4 gate is required\\n"); MPI_Abort(MPI_COMM(simpar),99); }\n    if(tx.on && postStage!=periodic_postStage_blend) MPI_Abort(MPI_COMM(simpar),99);\n    GfsWallEventStep wall;\n    gfs_wall_begin(simpar,&wall,postStage);',
        '    if(tx.on && postStage!=periodic_postStage_blend) MPI_Abort(MPI_COMM(simpar),99);\n    GfsWallEventStep wall;\n    gfs_wall_begin(simpar,&wall,postStage);\n    if(thermal64 && !tx.on && !wall.on){ fprintf(stderr,"[THERMAL64] transactional RK4 gate is required\\n"); MPI_Abort(MPI_COMM(simpar),99); }')
    replace_once(c, '    if(thermal64) gfs_t64_eos(simpar);\n    if(tx.on){',
        '    if(thermal64 && wall.on && !gfs_wall_valid(simpar,&wall,0)) goto gfs_wall_trial_rejected;\n    if(thermal64) gfs_t64_eos(simpar);\n    if(tx.on){')
    for name, cmd in [('exam', ['make', '-C', 'Exam', '-j2', 'USE_CUDA=1', 'GFS_VARPRESS_W0=1', 'libexam.a']),
                      ('link', ['make', '-j2', 'USE_CUDA=1', 'GFS_VARPRESS_W0=1', 'this'])]:
        with (root/('build_'+name+'.log')).open('x') as log:
            subprocess.run(cmd, cwd=dest, stdout=log, stderr=subprocess.STDOUT, check=True)
    result = {'host_sha256': sha(c), 'wall_sha256': sha(wall), 'thermal_sha256': sha(thermal),
              'binary_sha256': sha(dest/'eunha2.exe'),
              'scope': 'Experimental bounded one-rank W0 only. Build is not physics, rollback, wall impact, conservation, restart, or convergence validation.'}
    with (root/'build_result.json').open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps(result, indent=2))

if __name__ == '__main__': build(Path(sys.argv[1]).resolve())
