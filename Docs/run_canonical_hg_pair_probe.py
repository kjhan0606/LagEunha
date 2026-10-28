"""Dump t6 HG modes once with the existing verified single-rank probe binary."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import argparse
import struct

base=Path('/gpfs/kjhan/LagForce')
ap=argparse.ArgumentParser();ap.add_argument('--epoch',choices=['t0','t2','t6'],default='t6');args=ap.parse_args()
epoch=args.epoch
root=base/('Codex_pair_hg_modes_20261018' if epoch=='t6' else f'Codex_pair_hg_modes_{epoch}_20261018')
source={'t0':base/'Codex_canonical_initial_20261017/run/khout.000000.dat',
        't2':base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat',
        't6':base/'Codex_HG_constrained_20261017/omp8_t20/khout.000656.dat'}[epoch]
template=base/'Codex_mpi_hg_modes_20261017/np1'
binary=base/'Codex_mpi_hg_modes_20261017/eunha2.exe'
expected='cd730d83c1e404ce3418c9751dbd71a37d58ffdf853ca86b3bba87844a2a1f8a'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(binary)==expected
expected_source={'t0':'001e3ebbdfe5d006b9a3c3e510ca644ae347c2610f76d168a99b773c3c6cf862',
                 't2':'cafde902076beac9f8ccd0efc5bdec0a2c6b07e35c5fd46736fac427f56cc840',
                 't6':'a5101ac80d32ba07ee53cec39f5e0147e4e71ee997b86451502b62e103d00386'}[epoch]
assert sha(source)==expected_source
raw=source.read_bytes();_,n,size,_=struct.unpack_from('<iiii',raw)
source_time=struct.unpack_from('<d',raw,16+n*size)[0]
# icont=0 regenerates IC in kh.c. A byte-identical filename alias with index1
# lets the reader load the exact t0 checkpoint. Header time remains zero;
# no RK step or physical state change is performed by this diagnostic.
restart_index={'t0':1,'t2':184,'t6':656}[epoch]
source_copy=f'khout.{restart_index:06d}.dat'
root.mkdir(exist_ok=False)
(root/source_copy).write_bytes(raw)
(root/'params.dat').write_bytes((template/'params.dat').read_bytes())
env=json.loads((template/'environment.json').read_text())
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',HYDRO_TSTOP=str(source_time+.001),
           GFS_SHF_DUMP_FILE=str(root/'raw_faces.dat'))
assert env['GFS_MPI_HG_MODES_ONLY']=='1' and env['GFS_MPI_GEOMETRY_ONLY']=='1'
assert 'GFS_MPI_HG_APPLY_TAU' not in env
(root/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
(root/'manifest.json').write_text(json.dumps({'source':str(source),'source_sha256':sha(source),
    'binary':str(binary),'binary_sha256':expected,'purpose':'Read-only HG mode/RHS assembly at saved state.',
    'epoch':epoch,'source_time':source_time,'source_copy':source_copy,'restart_index':restart_index,
    'time_evolution':False,'ranks':1,'threads':1},indent=2)+'\n')
effective={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))}
effective.update(env)
with (root/'log').open('x') as log:
    p=subprocess.Popen(['mpirun','-np','1',str(binary),'params.dat',str(restart_index)],cwd=root,env=effective,stdout=log,stderr=subprocess.STDOUT)
    (root/'runner_pid').write_text(str(p.pid)+'\n')
    print(json.dumps({'root':str(root),'pid':p.pid}),flush=True)
    code=p.wait()
(root/'exit.code').write_text(str(code)+'\n')
assert code==0,code
log=(root/'log').read_text()
assert '[MPI_HG_MODES]' in log and '[ACCEPTED_STATE]' not in log
print(json.dumps({'exit':code,'modes_bytes':(root/'raw_faces.dat.hgm.rank0000').stat().st_size}),flush=True)
