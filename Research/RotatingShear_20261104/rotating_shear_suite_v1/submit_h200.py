from pathlib import Path
import subprocess,json,re,socket
root=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/rotating_shear_suite_v1')
assert socket.gethostname().split('.')[0]=='syntax'
rows=[]
for name in ['gresho_128x128','kh_128x256','kh_256x512']:
 d=root/name
 assert not (d/'submission.json').exists()
 with (d/'submission.lock').open('x') as f:f.write('Submission initiated once; inspect Slurm before any retry.\n')
 r=subprocess.run(['sbatch','--parsable',str(d/'job.sbatch')],text=True,capture_output=True)
 (d/'submission.stdout').write_text(r.stdout);(d/'submission.stderr').write_text(r.stderr)
 row={'case':name,'submit_exit':r.returncode,'job_id':r.stdout.strip().split(';')[0] if r.returncode==0 else None,'dependency':None}
 if r.returncode==0:assert re.fullmatch(r'\d+',row['job_id'])
 (d/'submission.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row);print(json.dumps(row),flush=True)
(root/'submissions.json').write_text(json.dumps(rows,indent=2)+'\n')
