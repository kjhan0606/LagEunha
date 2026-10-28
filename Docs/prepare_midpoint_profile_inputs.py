"""Validate and record the saved full-provider replay input, without replay."""
from pathlib import Path
import hashlib,json
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_compiled_midpoint_once_20261018';manifest=json.loads((root/'manifest.json').read_text());stage0=np.load(root/'stage0/inputs.npz');mid=np.load(root/'midpoint/inputs.npz');half=np.load(root/'half_update/candidate.npz');initial=np.load(root/'initial.npz')
def sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
for key in ('state','factors','mass','ids','frozen'):assert np.array_equal(mid[key],half[key])
seed=stage0['weight_rates'];N=len(mid['ids']);a0=initial['state'][:-4].reshape(-1,20);ah=mid['state'][:-4].reshape(-1,20);expected=a0[:,2]+5e-5*seed;expected[mid['frozen'].astype(bool)]=a0[mid['frozen'].astype(bool),2]
assert np.array_equal(ah[:,2],expected)
runner=base/'Codex_compiled_midpoint_modules_20261018/run_compiled_midpoint_once.py';assert hashlib.sha256(runner.read_bytes()).hexdigest()==manifest['runner_sha256']
modules={name:dict(path=meta['path'],sha256=hashlib.sha256(Path(meta['path']).read_bytes()).hexdigest()) for name,meta in manifest['modules'].items()};assert all(modules[k]['sha256']==meta['sha256'] for k,meta in manifest['modules'].items())
report=dict(status='prepared_not_executed',scope=__doc__,state_path=str(root/'midpoint/inputs.npz'),initial_seed_path=str(root/'stage0/inputs.npz'),initial_seed_key='weight_rates',time=float(mid['time']),seed_sha256=sha(seed),state_sha256=sha(mid['state']),factor_sha256=sha(mid['factors']),frozen_sha256=sha(mid['frozen']),half_checkpoint_equals_provider_input=True,half_weights_equal_initial_plus_half_dt_times_seed=True,runner=dict(path=str(runner),sha256=manifest['runner_sha256']),modules=modules,seed_contract='Original runner assigns d0,o0,g0,z0=evaluate(stage0), then evaluate(h/2,half,halfL,z0,midpoint). Thus stage0 converged weight_rates are the midpoint OUTER initial seed. The final midpoint weight_rates are not its initial seed.',outer=dict(maximum_iterations=10,termination_weight_rate_change=1e-11,required_mass_residual=2e-13),performance_replay_authorization='One complete provider replay only after fine completion, with same original inputs, flags, seed and stopping rule. No timestep update or tuning reruns. Instrumentation code not yet executed.',limits='Prepared inputs only, not a profiling or physics result.')
print(json.dumps(report,indent=2))
