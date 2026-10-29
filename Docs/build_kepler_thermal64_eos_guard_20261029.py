"""Build isolated EOS overflow/underflow trial guard for bounded thermal64."""
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
def replace_once(path,old,new):
    t=path.read_text()
    if t.count(old)!=1:raise ValueError('Unexpected context '+old[:100])
    path.write_text(t.replace(old,new));os.utime(path,None)
source=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallCandidate_20261029/project')
if sha(source/'Exam/exam.c')!='068a0db5638e908ea2b606bd2fdc7f4634263189b22185753af548263ef84f27' or sha(source/'Exam/gfs_wall_events.inc')!='7a831bfac2a6ac4c1823f253309f1f0d31835899ea76fbb2e95c139c03138be4':raise ValueError('Changed parent source')
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
wall=exam/'gfs_wall_events.inc';c=exam/'exam.c'
helper='''/* Reject an unrepresentable trial before gfs_t64_eos can abort.
 * No state, pressure floor or internal energy is modified here. */
static int gfs_wall_t64_eos_valid(SimParameters *s,GfsWallEventStep *w){
    if(!w->on || !gfs_thermal64_requested())return 1;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(s);
    for(ptrdiff_t i=0;i<w->np;i++){
        double U=GFS_T64_U(p[i]),V=GFS_T64_V(p[i]);
        double rho=p[i].mass/V,P=(GAS_GAMMA(s)-1.)*U/V;
        double cs=sqrt(GAS_GAMMA(s)*P/rho);
        float rf=(float)rho,pf=(float)P,cf=(float)cs;
        if(!GFS_T64_TAGGED(p[i]) || !(U>0) || !(V>0) || !isfinite(U) || !isfinite(V) ||
           !(rho>0) || !(P>0) || !(cs>0) || !isfinite(rho) || !isfinite(P) || !isfinite(cs) ||
           !(rf>0) || !(pf>0) || !(cf>0) || !isfinite(rf) || !isfinite(pf) || !isfinite(cf)){
            w->reason="nonpositive_or_nonfinite_trial";return 0;
        }
    }
    return 1;
}
'''
marker='static int gfs_wall_endpoint_inside(SimParameters *s,GfsWallEventStep *w,double h){'
replace_once(wall,marker,helper+marker)
old='    if(thermal64 && wall.on && !gfs_wall_valid(simpar,&wall,0)) goto gfs_wall_trial_rejected;\n    if(thermal64) gfs_t64_eos(simpar);'
new='    if(thermal64 && wall.on && (!gfs_wall_valid(simpar,&wall,0) || !gfs_wall_t64_eos_valid(simpar,&wall))) goto gfs_wall_trial_rejected;\n    if(thermal64) gfs_t64_eos(simpar);'
replace_once(c,old,new)
for name,cmd in [('exam',['make','-C','Exam','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','libexam.a']),('link',['make','-j2','USE_CUDA=1','GFS_VARPRESS_W0=1','this'])]:
    with (root/('build_'+name+'.log')).open('x') as f:subprocess.run(cmd,cwd=dest,stdout=f,stderr=subprocess.STDOUT,check=True)
result={'host_sha256':sha(c),'wall_sha256':sha(wall),'binary_sha256':sha(dest/'eunha2.exe'),
        'scope':'Experimental guard candidate only. Accepted-start initialization and other failure paths are unchanged; native rollback and physics parity require validation.'}
with (root/'build_result.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
