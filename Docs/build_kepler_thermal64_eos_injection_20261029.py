"""Build a verification-only EOS invalid-trial injection executable."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
source=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64EOSGuard_20261029/project')
if sha(source/'Exam/exam.c')!='e9117afaf36cff557063f802da3eb4d125dcb24ee1ca77e8712cf28edcbbaecc':raise ValueError('Changed guarded source')
root=Path(sys.argv[1]).resolve();dest=root/'project'
if dest.exists():raise FileExistsError(dest)
dest.mkdir(parents=True)
for p in source.iterdir():
    if p.name in ('Exam','eunha2.exe'):continue
    if p.name=='eunha2.o':shutil.copy2(p,dest/p.name)
    else:(dest/p.name).symlink_to(p.resolve(),target_is_directory=p.is_dir())
exam=dest/'Exam';exam.mkdir()
for p in (source/'Exam').iterdir():
    if p.is_dir():(exam/p.name).symlink_to(p.resolve(),target_is_directory=True)
    elif p.is_file():shutil.copy2(p,exam/p.name)
c=exam/'exam.c';t=c.read_text()
marker='    if(thermal64 && wall.on && (!gfs_wall_valid(simpar,&wall,0) || !gfs_wall_t64_eos_valid(simpar,&wall))) goto gfs_wall_trial_rejected;'
if t.count(marker)!=1:raise ValueError('Unexpected guard')
inject='''    /* Verification-only executable: corrupt one trial U once, never accepted. */
    static int eos_injection_done=0;
    if(thermal64 && wall.on && !eos_injection_done && getenv("GFS_TEST_EOS_INVALID_ONCE")){
        treevorostressrk4particletype *probe=(void*)VORORK4_TBP(simpar);
        GFS_T64_U(probe[0])=DBL_MAX;eos_injection_done=1;
        fprintf(stderr,"[EOS_INJECT_TEST] invalid trial injected before EOS guard\\n");
    }
'''
c.write_text(t.replace(marker,inject+marker));os.utime(c,None)
for name,cmd in [('exam',['make','-C','Exam','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','libexam.a']),('link',['make','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','this'])]:
    with (root/('build_'+name+'.log')).open('x') as f:subprocess.run(cmd,cwd=dest,stdout=f,stderr=subprocess.STDOUT,check=True)
result={'host_sha256':sha(c),'binary_sha256':sha(dest/'eunha2.exe'),
        'scope':'Verification-only executable containing deliberate one-time test corruption; must never be used for production.'}
with (root/'build_result.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
