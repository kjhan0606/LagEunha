from pathlib import Path
import json, os, subprocess, time
rank_text=os.environ.get('NATIVE_MPI_RANKS') or os.environ.get('SLURM_NTASKS')
if rank_text is None:
    raise RuntimeError('Set NATIVE_MPI_RANKS explicitly for a manual native launch')
ranks=int(rank_text)
if ranks<1:
    raise ValueError('Native MPI rank count must be positive')
t0=time.monotonic()
with open('log','w') as out:
    p=subprocess.run(['/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun','-np',str(ranks),'./eunha2.exe','params.dat'] + ([os.environ['KEPLER_RESTART_STEP']] if 'KEPLER_RESTART_STEP' in os.environ else []),stdout=out,stderr=subprocess.STDOUT)
    out.write(f'EXIT:{p.returncode}\n')
Path('timing.json').write_text(json.dumps({'elapsed_seconds':time.monotonic()-t0,'ranks':ranks,'exit_code':p.returncode},indent=2)+'\n')
raise SystemExit(p.returncode)
