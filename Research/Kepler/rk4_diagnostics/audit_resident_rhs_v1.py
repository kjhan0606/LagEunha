"""Independent saved-state/mass/rotation audit of the actual next RK k1."""
import csv
import hashlib
import json
import pathlib
import struct
from decimal import Decimal, getcontext
import numpy as np

getcontext().prec=60
root=pathlib.Path(__file__).resolve().parent
source=root/'material_full_resident_rhs_v1_ada_mpi4'
dest=root/'material_full_resident_rhs_audit_v1'
dest.mkdir()
assert (source/'native_exit_code.txt').read_text().strip()=='0'
next_rhs=json.loads((source/'run/next_RHS_result.json').read_text())
rhsdir=source/'run/next_RK_stage1'
rhs=json.loads((rhsdir/'rhs.json').read_text(),parse_float=Decimal)
mass=json.loads((rhsdir/'GPU_mass.json').read_text(),parse_float=Decimal)
resident=json.loads((rhsdir/'resident_mass.json').read_text(),parse_float=Decimal)
assert rhs['RHS_gates_pass'] and max(mass['true_relative_residual'])<Decimal('1e-16')
assert resident['same_RHS_old_host_vector_DD_reference']
assert resident['reference_work_difference_over_U']<Decimal('1e-8')
assert resident['reference_torque_difference_over_2K']<Decimal('1e-11')
assert next_rhs['accepted_steps_unchanged']==4 and not next_rhs['new_accepted_RK_step']
checked=0
for name in ('source.sha256','binary.sha256','accepted_input.sha256'):
 for line in (source/name).read_text().splitlines():
  digest,path=line.split(maxsplit=1)
  assert hashlib.file_digest(pathlib.Path(path).open('rb'),'sha256').hexdigest()==digest
  checked+=1
seen=set();energy=angular=Decimal(0)
for rank in range(4):
 with (rhsdir/f'cells_r{rank}.csv').open() as f:
  for row in csv.DictReader(f):
   cell=int(row['cell']);assert cell not in seen;seen.add(cell)
   assert Decimal(row['volume'])>0 and Decimal(row['density'])>0 and Decimal(row['U'])>0
   energy=max(energy,Decimal(row['energy_rate_error_over_U']))
   angular=max(angular,Decimal(row['angular_rate_error_over_2K']))
   spin,inertia,rot=map(Decimal,(row['spin'],row['inertia'],row['Erot']))
   assert inertia>0 and abs(rot-spin*spin/(2*inertia))<Decimal('1e-18')*max(rot,Decimal('1e-40'))
 with (rhsdir/f'prepared_k1_r{rank}.bin').open('rb') as f:
  h=np.fromfile(f,'<u8',8);assert int(h[0])==0x4d464e4558544b31 and int(h[1])==1 and int(h[2])==rank
  time=struct.unpack('<d',f.read(8))[0];assert time==next_rhs['accepted_time_unchanged']
  dt=np.fromfile(f,np.longdouble,1);rates=np.fromfile(f,np.longdouble,3)
  k=np.fromfile(f,'<f8',4*int(h[4])).reshape(-1,4);assert not f.read() and np.isfinite(k).all()
 checkpoint=root/'material_full_cached_continuation_v1_ada_mpi4/run/continuous_step4'/f'accepted_r{rank}.bin'
 with checkpoint.open('rb') as f:
  old=np.fromfile(f,'<u8',8);assert int(old[1])==2 and old[4]==h[4] and old[5]==32
  length=struct.unpack('<Q',f.read(8))[0];f.seek(length+int(old[6])+10*int(old[7]),1)
  accepted=np.fromfile(f,'<f8',4*int(h[4])).reshape(-1,4)
 assert np.array_equal(k[:,:2].view('<u8'),accepted[:,2:].view('<u8'))
assert len(seen)==71080 and energy<Decimal('1e-8') and angular<Decimal('1e-11')
out={'scope':'saved actual resident-mass k1, all 71080 source ledgers and exact saved-state xdot=v check',
     'audit_pass':True,'hashes_checked':checked,'source_cells_checked':len(seen),
     'energy_rate_error_over_U':str(energy),'angular_rate_error_over_2K':str(angular),
     'primary_solve_seconds':str(resident['primary_solve_seconds_rank0']),
     'same_RHS_reference_seconds':str(resident['reference_solve_seconds_rank0']),
     'prepared_k1_matches_accepted_velocity_bits':True,'accepted_time_unchanged':next_rhs['accepted_time_unchanged'],
     'new_accepted_RK_step':False,'long_Kepler_validated':False}
(dest/'contract.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
