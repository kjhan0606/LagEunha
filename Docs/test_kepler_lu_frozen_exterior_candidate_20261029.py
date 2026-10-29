import importlib.util, os
from pathlib import Path
import numpy as np
from test_frozen_moment_boundary import state
import kinetic_factor_rhs as candidate
source=Path('/gpfs/kjhan/LagForce/Codex_LU_kepler_preflight_n64_20261017_093046/kinetic_factor_rhs.py')
spec=importlib.util.spec_from_file_location('old_kinetic_factor_rhs',source)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
for backend in ('scalar','batched'):
 os.environ['LAGFORCE_VECTOR_MOMENT_FLUX']='1' if backend=='batched' else '0'
 m,y,L=state();pos=y[:-4].reshape(-1,20)[:,0]
 for label,mask in [('none',np.zeros(len(m),bool)),('mixed',pos<.3),('all',np.ones(len(m),bool))]:
  a,da=old.rhs(0.,y,m,covariance_factors=L,frozen_mask=mask)
  b,db=candidate.rhs(0.,y,m,covariance_factors=L,frozen_mask=mask)
  delta=float(np.max(abs(a-b)));active=~mask
  areadiff=float(np.max(abs(da['area_rates'][active]-db['area_rates'][active]),initial=0.))
  if delta>2e-13 or areadiff>2e-13:raise AssertionError((backend,label,delta,areadiff))
  print(backend,label,'max_rhs_difference',delta,'active_area_rate_difference',areadiff,flush=True)
