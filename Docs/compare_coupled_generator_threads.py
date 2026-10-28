"""Compare saved OMP1/OMP4 one-step results without rerunning evolution."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');one=b/'Codex_coupled_generator_omp1_20261018';four=b/'Codex_coupled_generator_steps_20261018'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--one',type=Path,default=one);parser.add_argument('--four',type=Path,default=four)
parser.add_argument('--first-label',default='omp1');parser.add_argument('--second-label',default='omp4')
parser.add_argument('--output',default='thread_comparison.json')
args=parser.parse_args();one=args.one;four=args.four
out={}
for sign in ['plus','minus']:
 name=sign+'_5e-05'
 for root in [one,four]:assert (root/name/'exit.code').read_text().strip()=='0'
 p=one/name/'khout.000744.dat';q=four/name/'khout.000744.dat'
 a,t=read(p);c,s=read(q);assert t==s and a.dtype==c.dtype
 out[sign]={'time':t,'fields_bitwise_equal':{k:bool(np.array_equal(a[k],c[k])) for k in a.dtype.names},args.first_label+'_sha256':sha(p),args.second_label+'_sha256':sha(q)}
out['scope']=f'{args.first_label} versus {args.second_label}: two finite-perturbation branches, one accepted dt5e-5 step at t3.4. Reader fields compared; no long MPI or full CUDA-hydro claim.'
with (one/args.output).open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
