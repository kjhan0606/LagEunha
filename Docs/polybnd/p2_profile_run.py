"""Short isolated performance/physics comparison; never import running-job state."""
import os,sys,time,runpy,cProfile,pstats,json
from pathlib import Path
import numpy as np
mode=sys.argv[1];n=int(sys.argv[2]);end=float(sys.argv[3]);profile='--profile' in sys.argv
source='stagP2_perf_reference.py' if mode=='reference' else 'stagP2_perf_candidate.py'
os.environ['P2_GEOMETRY']=mode
sys.argv=[source,str(n),'.01',str(end),'2.78e-5','.84','23','0.0',mode,'.01']
p=cProfile.Profile();start=time.perf_counter()
if profile:p.enable()
g=runpy.run_path(source,run_name='__main__')
if profile:p.disable();p.dump_stats(mode+'.pstats');pstats.Stats(p).sort_stats('cumulative').print_stats(20)
elapsed=time.perf_counter()-start
if abs(g['t']-end)>1e-10:raise RuntimeError('short run did not reach target')
keys=sorted(g['zstate']);values=np.array([g['zstate'][k] for k in keys])
np.savez(mode+'.npz',X=g['X'],V=g['Vv'],keys=np.array(keys),zstate=values,mA=g['mA'],t=g['t'],ledgers=np.array([g[k] for k in ['S_tot','Eremap','nremap_tot','Ab_tot','nrej']]))
threads=None
if mode=='openmp':
    from p2_geometry_backend import threads_used
    threads=threads_used()
Path(mode+'.timing.json').write_text(json.dumps(dict(mode=mode,n=n,end=end,seconds=elapsed,profile=profile,omp_requested=os.environ.get('OMP_NUM_THREADS'),omp_actual=threads),indent=2))
