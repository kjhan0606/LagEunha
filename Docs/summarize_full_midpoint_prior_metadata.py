"""Numerical prior cost diagnostics on saved internal predictors, no RHS."""
import argparse,json,time
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
base=Path('/gpfs/kjhan/LagForce');reports={};contrasts={};densities={};started=time.perf_counter()
def quantile(x):return np.quantile(x,[0,.25,.5,.75,.9,.99,1]).tolist() if len(x) else []
for label,dirname in [('initial','Codex_full_disk_initial_density_20261019'),('coarse','Codex_full_disk_coarse_midpoint_density_20261019'),('fine','Codex_full_disk_fine_midpoint_density_20261019')]:
 rows=[json.loads(line) for line in (base/dirname/'cells.jsonl').read_text().splitlines()];r=np.load(base/dirname/'targets.npz')['raw'];assert all(row['status']=='passed' for row in rows)
 contrast=np.array([2*(x['prior']['high']-x['prior']['low'])/(x['prior']['high']+x['prior']['low']) if not x['prior']['flat'] else 0. for x in rows]);density=r[:,0]/np.array([x['parameter']['area'] for x in rows]);contrasts[label]=contrast;densities[label]=density
 reports[label]=dict(nonflat=int((contrast>0).sum()),relative_contrast_quantiles=quantile(contrast[contrast>0]))
mask=(contrasts['coarse']>0)&(contrasts['fine']>0);ratio=contrasts['coarse'][mask]/contrasts['fine'][mask]
change={k:(densities[k]-densities['initial'])/densities['initial'] for k in ('coarse','fine')}
result=dict(scope=__doc__,reports=reports,common_nonflat=int(mask.sum()),coarse_fine_contrast_ratio_quantiles=quantile(ratio),relative_cell_mean_changes={k:dict(quantiles=quantile(v),rms=float(np.sqrt(np.mean(v*v)))) for k,v in change.items()},provider_calls=0,evolution_steps=0,seconds=time.perf_counter()-started,limits='Internal forward predictors and numerical priors only; these ratios are not an endpoint convergence order or physical density-error measurement. No prior threshold or reconstruction changed.')
(args.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
