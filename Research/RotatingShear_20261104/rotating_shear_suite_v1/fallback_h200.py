from pathlib import Path
import subprocess,json,re
root=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/rotating_shear_suite_v1');rows=[]
for name in ['gresho_128x128','kh_128x256','kh_256x512']:
 d=root/name;meta=json.loads((d/'submission.json').read_text());job=meta['job_id']
 q=subprocess.run(['scontrol','show','job','-o',job],text=True,capture_output=True,check=True);(d/'fallback_before_cancel.txt').write_text(q.stdout)
 assert 'UserId=kjhan(' in q.stdout and f'WorkDir={d}' in q.stdout
 assert 'JobState=PENDING' in q.stdout and not (d/'process.json').exists()
 subprocess.run(['scancel','--state=PENDING',job],check=True)
 q=subprocess.run(['scontrol','show','job','-o',job],text=True,capture_output=True,check=True);(d/'fallback_after_cancel.txt').write_text(q.stdout)
 state=re.search(r'JobState=(\S+)',q.stdout).group(1);row={'case':name,'job_id':job,'state':state,'reason':'QOSGrpGRES blocks H200; withdraw only our never-started submissions and use authorized Ada'}
 assert state=='CANCELLED',row
 (d/'fallback.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row);print(json.dumps(row),flush=True)
(root/'fallbacks.json').write_text(json.dumps(rows,indent=2)+'\n')
