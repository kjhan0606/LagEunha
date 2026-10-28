"""Compare extracted cell preparation to unchanged whole-state old oracle."""
import argparse,hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_full_disk_density_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_face_state_batch_api_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018')]
os.environ.update(LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
from quadratic_cell_prepare import prepare_cell
import quadratic_velocity_reconstruction as quad
from saved_contact_density_restore import restore_density
import moment_ale_evolve as model
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
archive=base/'Codex_numerical_Psi_endpoints_v2_20261018';data=np.load(archive/'initial_reconstruction.npz');initial=np.load(base/'Codex_compiled_midpoint_once_20261018/initial.npz');mass=data['mass'];L=initial['factors'];ids=data['ids'];raw=np.column_stack([mass,data['state'][:-4].reshape(-1,20)[:,3:]])
priors=json.loads((archive/'initial_priors.json').read_text());fits=json.loads((archive/'initial_fits.json').read_text());off=data['polygon_offsets'];vertices=data['polygon_vertices'];polys=[vertices[a:b] for a,b in zip(off[:-1],off[1:])]
model.GAMMA=float(np.float32(1.666667));densities=[restore_density(p,r,prior,fit) for p,r,prior,fit in zip(polys,raw,priors,fits)];states=[model.decode(p,r,density=d) for p,r,d in zip(polys,raw,densities)]
for s,factor in zip(states,L):s['R']=factor@factor.T
old=quad.prepare(states,polys,[],mass,L);centers=np.array([s['c'] for s in states]);means=np.array([s['v'] for s in states]);second=np.asarray([s['raw'][3:7].reshape(2,2)/m for s,m in zip(states,mass)]);tree=cKDTree(centers);maximum=0.
for i,s in enumerate(states):
 metadata,diag=prepare_cell(i,s,densities[i],mass[i],L[i],centers,means,second,tree)
 for name,reference in s['velocity_quadratic'].items():
  got=np.asarray(metadata[name]);want=np.asarray(reference);maximum=max(maximum,float(abs(got-want).max(initial=0.)));assert np.array_equal(got,want),(int(ids[i]),name)
report=dict(status='passed',cells=len(ids),metadata_bitwise_equal=True,maximum_difference=maximum,original_oracle_diagnostics=old,seconds=time.perf_counter()-started,provider_calls=0,evolution_steps=0,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'quadratic_cell_prepare.py',Path(quad.__file__)]},limits='Stored initial1590 state, old fitted densities restored without Newton; no actual full-domain provider or evolution.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
