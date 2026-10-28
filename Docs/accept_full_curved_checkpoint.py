"""Preserve a validated numerical step as an exact restart checkpoint.

Acceptance here denotes factor/geometry/density gates, not long-time physics.
"""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
for k in ('state','factor','geometry','density','output'):parser.add_argument('--'+k,type=Path,required=True)
args=parser.parse_args();sha=hashlib.sha256(args.state.read_bytes()).hexdigest();saved=np.load(args.state);initial=np.load('/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/initial.npz');N=len(saved['ids']);reports={k:json.loads((getattr(args,k)/'report.json').read_text()) for k in ('factor','geometry','density')}
assert reports['factor']['status']=='passed_full_disk_factor_candidate'
assert reports['geometry']['status']=='passed_full_curved_stage_geometry' and reports['geometry']['state_sha256']==sha
assert reports['density']['status']=='passed_full_stage_Psi_density_preparation' and reports['density']['state_sha256']==sha and reports['density']['counts']=={'passed':N}
assert reports['density']['time']==reports['geometry']['time']==float(saved['time'])
assert np.array_equal(saved['ids'],initial['ids']) and np.array_equal(saved['mass'],initial['mass']) and not saved['frozen'].any() and saved['psi_constrained_mask'].all()
assert np.isfinite(saved['state']).all() and np.isfinite(saved['factors']).all() and np.isfinite(saved['Psi']).all()
args.output.mkdir(exist_ok=False);target=args.output/'accepted.npz';shutil.copy2(args.state,target);assert hashlib.sha256(target.read_bytes()).hexdigest()==sha
manifest=dict(status='accepted_numerical_checkpoint',scope=__doc__,time=float(saved['time']),cells=N,checkpoint=str(target),checkpoint_sha256=sha,source_checkpoint=str(args.state),source_reports={k:str(getattr(args,k)/'report.json') for k in reports},provider_calls=0,new_evolution_steps=0,copy_policy='Byte-identical checkpoint preservation; no reset/refit/projection/energy or mass correction',physical_assessment='Actual-cell reference errors and common-time dt comparison are reported separately. Long-term Kepler stability remains unproven.')
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
