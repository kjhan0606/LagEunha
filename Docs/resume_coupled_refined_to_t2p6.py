"""Continue completed amplitude/dt checks, preserving intermediate dumps."""
import argparse,json,os,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,default=b/'Codex_coupled_convergence_20261018')
parser.add_argument('--root',type=Path,default=b/'Codex_coupled_refined_to_t2p6_20261018')
parser.add_argument('--coarse-step',type=int,default=423)
args=parser.parse_args();src=args.source;root=args.root
assert (src/'analysis.json').exists(), 'Require completed source convergence analysis before continuation'
assert not root.exists(), 'Preserve existing results'
old=json.loads((src/'manifest.json').read_text());_,target=read(b/f'Codex_coupled_onset_to_t4_20261018/reference/khout.{args.coarse_step:06d}.dat')
# Validate every input before launching any branch. Source indices may
# already include an earlier continuation; never assume all start at184.
prepared={}
for name,info in old['cases'].items():
 assert (src/name/'exit.code').read_text().strip()=='0'
 step=info.get('source_step',184)+info['expected_steps']
 snap=f'khout.{step:06d}.dat';a,t=read(src/name/snap)
 steps=round((target-t)/info['dt'])
 assert steps>0 and abs(t+steps*info['dt']-target)<1e-10
 c=json.loads((src/name/'environment.json').read_text())
 assert (src/name/'params.dat').is_file()
 prepared[name]=(step,snap,t,steps,c)
dump_count=int((target-min(v[2] for v in prepared.values()))/.2)+2
estimated_bytes=sum((src/name/item[1]).stat().st_size*(dump_count+1) for name,item in prepared.items())
free_bytes=shutil.disk_usage(root.parent).free
assert free_bytes>2*estimated_bytes
print(json.dumps({'source':str(src),'destination':str(root),'target':target,
 'dump_interval':.2,'final_dump':True,'estimated_snapshot_bytes':estimated_bytes,
 'free_bytes':free_bytes}),flush=True)
root.mkdir(exist_ok=False)
meta={'target':target,'source_root':str(src),'cases':{}};jobs=[]
for name,info in old['cases'].items():
 step,snap,t,steps,c=prepared[name]
 d=root/name;d.mkdir()
 for f in ['params.dat',snap]:shutil.copyfile(src/name/f,d/f)
 c.update(HYDRO_TSTOP=str(target-.5*min(v['dt'] for v in old['cases'].values())),EUNHA_DUMP_DT='.2',EUNHA_DUMP_FINAL='1')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','5400','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat',str(step)],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);jobs.append((name,d,p,f))
 meta['cases'][name]={**info,'pid':p.pid,'source_step':step,'start_time':t,'expected_steps':steps};(root/'manifest.json').write_text(json.dumps(meta,indent=2));print(name,p.pid,flush=True)
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(name,code,flush=True)
