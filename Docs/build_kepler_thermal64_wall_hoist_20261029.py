"""Build experimental wall candidate hoisting stable thermal flag per traversal."""
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
def change(text,old,new):
    if text.count(old)!=1:raise ValueError('Unexpected wall source: '+old[:80])
    return text.replace(old,new)
source=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallCandidate_20261029/project')
if sha(source/'Exam/gfs_wall_events.inc')!='7a831bfac2a6ac4c1823f253309f1f0d31835899ea76fbb2e95c139c03138be4':raise ValueError('Changed base wall source')
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
wall=exam/'gfs_wall_events.inc';text=wall.read_text()
text=change(text,'long double a[7]={0,0,0,0,0,0,0};\n    for(ptrdiff_t i=0;i<VORO_NP(s);i++){',
    'long double a[7]={0,0,0,0,0,0,0};\n    const int thermal64=gfs_thermal64_requested();\n    for(ptrdiff_t i=0;i<VORO_NP(s);i++){')
text=change(text,'double U=gfs_t64_U(p+i);','double U=thermal64 && GFS_T64_TAGGED(p[i])?GFS_T64_U(p[i]):(double)p[i].ie;')
text=change(text,'    for(ptrdiff_t i=0;i<w->np;i++){\n        size_t id=(size_t)PINDX((treevorork4particletype*)(p+i));\n        if((gfs_thermal64_requested()',
    '    const int thermal64=gfs_thermal64_requested();\n    for(ptrdiff_t i=0;i<w->np;i++){\n        size_t id=(size_t)PINDX((treevorork4particletype*)(p+i));\n        if((thermal64')
text=change(text,'    gfs_wall_new_generation(w);\n    for(ptrdiff_t i=0;i<w->np;i++){\n        /* N distinct members',
    '    gfs_wall_new_generation(w);\n    const int thermal64=gfs_thermal64_requested();\n    for(ptrdiff_t i=0;i<w->np;i++){\n        /* N distinct members')
text=change(text,'if((gfs_thermal64_requested() && (!GFS_T64_TAGGED(p[i]) || !isfinite(GFS_T64_U(p[i]))',
    'if((thermal64 && (!GFS_T64_TAGGED(p[i]) || !isfinite(GFS_T64_U(p[i]))')
wall.write_text(text);os.utime(wall,None);os.utime(exam/'exam.c',None)
for name,cmd in [('exam',['make','-C','Exam','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','libexam.a']),('link',['make','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','this'])]:
    with (root/('build_'+name+'.log')).open('x') as log:subprocess.run(cmd,cwd=dest,stdout=log,stderr=subprocess.STDOUT,check=True)
result={'wall_source_sha256':sha(wall),'host_source_sha256':sha(exam/'exam.c'),'binary_sha256':sha(dest/'eunha2.exe'),
        'scope':'Experimental flag hoist in bounded wall traversals only; no physics, parity or performance qualification implied.'}
with (root/'build_result.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
