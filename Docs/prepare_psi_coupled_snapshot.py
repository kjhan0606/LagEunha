"""Reproduce a fresh immutable two-provider connection snapshot.

Run editing/small-file preparation only; refuses existing destination.
This does not execute any calculation. Tested historical snapshots stay intact.
"""
import argparse,hashlib,json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--parent',type=Path,default=Path('/gpfs/kjhan/LagForce/Codex_face_record_native_modules_v2_20261019'))
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();repo=Path(__file__).resolve().parent
parent=json.loads((args.parent/'manifest.json').read_text())
# Verify every inherited source before making a runnable copy.
for name,digest in parent['python_modules'].items():assert hashlib.sha256((args.parent/(name+'.py')).read_bytes()).hexdigest()==digest
for name,digest in parent['native_files'].items():assert hashlib.sha256((args.parent/name).read_bytes()).hexdigest()==digest
args.output.mkdir(exist_ok=False)
for name in parent['python_modules']:(args.output/(name+'.py')).write_bytes((args.parent/(name+'.py')).read_bytes())
for name in parent['native_files']:(args.output/name).write_bytes((args.parent/name).read_bytes())
for name in ('contact_potential_density.py','positive_radial_quadrature.py','circular_power_reference.py'):(args.output/name).write_bytes((repo/'polybnd'/name).read_bytes())
path=args.output/'kinetic_factor_rhs.py';source=path.read_text();a=source.index('    if audit_weight_rates is not None:',source.index('def rhs('));b=source.index('    dr=np.zeros((N,18));ledger=',a)
source=source[:a]+'    coupled_root_begin=time.perf_counter();coupled_root_cpu=time.process_time()\n'+source[a:b]+'    coupled_massroot_timing=dict(wall_seconds=time.perf_counter()-coupled_root_begin,process_cpu_seconds=time.process_time()-coupled_root_cpu)\n'+source[b:]
source=source.replace('    diag.update(face_record_total_timing=record_timing,',"    diag['coupled_massroot_timing']=coupled_massroot_timing\n    diag.update(face_record_total_timing=record_timing,")
path.write_text(source)
manifest=dict(parent_snapshot=str(args.parent),parent_manifest_sha256=hashlib.sha256((args.parent/'manifest.json').read_bytes()).hexdigest(),python_modules={},native_files={},scope='Psi-coupled API plus one large massroot timer; inherited arithmetic unchanged',preparer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
for name in list(parent['python_modules'])+['contact_potential_density','positive_radial_quadrature','circular_power_reference']:manifest['python_modules'][name]=hashlib.sha256((args.output/(name+'.py')).read_bytes()).hexdigest()
for name in parent['native_files']:manifest['native_files'][name]=hashlib.sha256((args.output/name).read_bytes()).hexdigest()
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(dict(status='prepared_only',snapshot=str(args.output))))
