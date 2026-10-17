"""One-shot Slurm/checkpoint monitor; caller must space invocations >=60s."""
import json
import subprocess
from pathlib import Path
from datetime import datetime,timezone

RUN=Path('/gpfs/kjhan/LagForce/Codex_P_adaptive_coupled_H200_20261010_200459')
CONTROL=Path('/gpfs/kjhan/LagForce/Codex_P_parallel_20261010_222051')
WAVE=Path('/gpfs/kjhan/LagForce/Codex_P_wave_parallel_20261010_222051')

def snapshot():
    account=subprocess.run(['sacct','-j','418552,418553,418555','-X','--noheader','--parsable2',
        '--format=JobID,State,Elapsed,ExitCode'],capture_output=True,text=True,check=True)
    jobs=[]
    for row in account.stdout.splitlines():
        parts=row.strip().split('|')
        if len(parts)>=4:jobs.append(dict(zip(('job','state','elapsed','exit'),parts[:4])))
    cases={}
    for case in ('shear','wave'):
        p=(RUN if case=='shear' else WAVE)/case;entry={}
        for name in ('checkpoint.json','audit_checkpoint.json','result.json'):
            f=p/name
            if not f.exists():continue
            data=json.loads(f.read_text())
            entry[name]={k:v for k,v in data.items() if k not in
                         ('config','state','factors','mass','signed','absolute','quad_delta')}
        cases[case]=entry
    result=dict(checked_at=datetime.now(timezone.utc).isoformat(),jobs=jobs,cases=cases)
    tmp=CONTROL/'latest_status.tmp';tmp.write_text(json.dumps(result,indent=2));tmp.replace(CONTROL/'latest_status.json')
    print(json.dumps(result))

if __name__=='__main__':snapshot()
