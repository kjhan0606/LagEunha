"""Short MPI rank-parity checks for the already approved SSPRK2 implementation.
No new physics, IC model, Slurm jobs, or production binaries are introduced.
"""
from pathlib import Path
import subprocess,sys
O=Path(__file__).resolve().parent;P=O/'qualification.py';B=O/'build_final/native.exe'
tasks=[]
for case,source in [('sod','sod_step_dt00005'),('shear','shear_interior32_t01'),('gresho','gresho_interior64_t03')]:
 for rank in [1,4]:tasks.append([str(O.parent/source),case+'_final_r'+str(rank),'--ranks',str(rank),'--window','.000025','--step','1'])
for rank in [1,4]:tasks.append([str(O/'noh_r1_maps'),'noh_stages_r'+str(rank),'--ranks',str(rank),'--volume','stages','--window','.0001'])
for rank in [1,4]:tasks.append([str(O/'kh_r1_suite'),'kh_reject_verified_r'+str(rank),'--ranks',str(rank),'--binary',str(O/'build_reject/native.exe')])
for rank in [1,4]:tasks.append([str(O/'rt_r1_suite'),'rt_final_r'+str(rank),'--ranks',str(rank),'--window','.0001'])
for rank in [1,4]:tasks.append([str(O/'sedov_r1_suite'),'sedov_final_r'+str(rank),'--ranks',str(rank),'--window','.0001'])
for rank in [1,2,4]:tasks.append([str(O/'cyl_event_r4'),'cyl_restart_r'+str(rank),'--ranks',str(rank),'--window','.0001'])
for task in tasks:
 cmd=[sys.executable,str(P)]+task
 if '--binary' not in task:cmd+=['--binary',str(B)]
 subprocess.run(cmd,cwd=O,check=True)
print('FINAL_SUITE_COMPLETED',flush=True)
