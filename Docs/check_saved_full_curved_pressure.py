"""Preliminary candidate pressure on saved actual cell volumes; no density fit."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
for name in ('state','geometry','output'):parser.add_argument('--'+name,type=Path,required=True)
args=parser.parse_args();args.output.mkdir(exist_ok=False);base=Path('/gpfs/kjhan/LagForce');d=np.load(args.state);a=d['state'][:-4].reshape(-1,20);mesh=np.load(args.geometry/'mesh.npz');sha=hashlib.sha256(args.state.read_bytes()).hexdigest();assert json.loads((args.geometry/'report.json').read_text())['state_sha256']==sha
p=json.loads((base/'Codex_full_curved_kepler_IC_20261019/report.json').read_text())['physics'];thermal=a[:,19]-.5*(a[:,15]+a[:,18]);pressure=(p['gamma']-1)*thermal/mesh['moments'][:,0];error=(pressure-p['P0'])/p['P0'];initial=np.load(base/'Codex_full_curved_kepler_IC_20261019/initial.npz');r=np.linalg.norm(initial['state'][:-4].reshape(-1,20)[:,:2],axis=1);regions={}
for name,mask in dict(all=np.ones(len(r),bool),inner_initial_radius_lt_0p4=r<.4,disk_initial_radius_0p5_to_2=(r>.5)&(r<2)).items():
 e=error[mask];m=d['mass'][mask];regions[name]=dict(cells=int(mask.sum()),relative_pressure_RMS=float(np.sqrt(np.average(e*e,weights=m))),maximum_absolute_relative_error=float(abs(e).max()),minimum_pressure=float(pressure[mask].min()),maximum_pressure=float(pressure[mask].max()))
report=dict(status='completed_preliminary_candidate_pressure',scope=__doc__,state_path=str(args.state),state_sha256=sha,time=float(d['time']),regions=regions,provider_calls=0,evolution_steps=0,limits='Stored thermal energy divided by actual cell area against prescribed constant P0. Density realizability/full actual-cell field reference/two-dt and long-term checks remain separate. No reset or pressure correction.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
