"""Canonical ledger ABI against the original full-stage backend, same host."""
import argparse,json,sys,time,hashlib
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
sys.path.insert(0,str(args.library.parent));sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_compiled_midpoint_modules_20261018')
from moment_owned_backend import OwnedBackend
base=Path('/gpfs/kjhan/LagForce');c=base/'Codex_kinetic_face_capture_20261018';s=base/'Codex_moment_stage_capture_20261018';h=base/'Codex_hydro_cell_backends_20261018'
def read(root,name,dtype,width=None):
 a=np.fromfile(root/name,dtype=dtype);return a.reshape(-1,width) if width else a
data=dict(node_inputs=read(c,'inputs.bin','<f8',33),cells=read(h,'cells.bin','<f8',19),owners=read(h,'owners.bin','<i8',2),offsets=read(h,'offsets.bin','<i8'),incidence=read(h,'incidence.bin','<i8'),frozen=read(h,'frozen.bin','<i8'),extra=read(s,'node_extra.bin','<f8',20),factors=read(s,'factors.bin','<f8',4).reshape(-1,2,2),gravity_nodes=read(s,'gravity_nodes.bin','<f8',7),gravity_offsets=read(s,'gravity_offsets.bin','<i8'))
args.output.mkdir(exist_ok=False);rows=[];start=time.perf_counter()
for kind,threads in (('cpu',1),('cpu',4),('cuda',1)):
 backend=OwnedBackend(args.library,kind,threads);nf=backend.nodes(data['node_inputs']);out=backend.stage(data)
 led=backend.ordered_ledgers(data,out,nf);assert np.array_equal(led,out['ledgers'])
 failures={}
 for name in ('bad_status','bad_owner','nonfinite_flux'):
  d=dict(data);o=dict(out);f=nf
  if name=='bad_status':o['stage']=out['stage'].copy();o['stage'][len(data['cells'])//2,93]=1.
  elif name=='bad_owner':d['owners']=data['owners'].copy();d['owners'][0,0]=-1
  else:f=nf.copy();f[-1,0]=np.nan
  status,sentinel=backend.ordered_ledgers(d,o,f,sentinel=-947.)
  assert status and np.all(sentinel==-947.),(kind,name,status)
  failures[name]=dict(status=status,caller_output_unchanged=True)
 rows.append(dict(backend=kind,threads=threads,ledger_bitwise=True,failures=failures))
report=dict(status='passed_owned_ordered_ledger_component',cells=len(data['cells']),nodes=len(nf),backends=rows,library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),seconds=time.perf_counter()-start,physical_providers=0,evolution_steps=0,limits='Canonical ledger component only; no MPI evolution claim')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
