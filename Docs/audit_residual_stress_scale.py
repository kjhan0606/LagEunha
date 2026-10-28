"""Compare saved residual velocity covariance scale to thermal pressure; no closure adopted."""
import ast,json,importlib.util
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_hybrid_affine_path_20261018'
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
rd=load('rd',base/'run_canonical_av_branches_20261018.py');a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
tree=ast.parse((base/'audit_remap_true_spatial.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='recover'],type_ignores=[]),'recover','exec'),globals())
# Explicit comparison gamma, not a readback of current production parameters.
gamma=5/3;rows=[]
for N in (4,8):
 d=np.load(root/f'steps{N}.npz');q=d['state'];ids=d['ids'];S,C,G,R=recover(q);U=q[:,15]-.5*(q[:,12]+q[:,14]);ev=np.linalg.eigvalsh(R);r=np.hypot(a['x'][ids]-3,a['y'][ids]-3);regions={}
 for name,mask in [('inner',r<.5),('edge',(r>=.5)&(r<.6)),('patch',r<.7)]:
  ratio=ev[mask,1]/((gamma-1)*U[mask]);trace=np.trace(R[mask],axis1=1,axis2=2);positive=trace>1e-14*U[mask];anis=(ev[mask,1]-ev[mask,0])[positive]/trace[positive]
  regions[name]={'cells':int(mask.sum()),'residual_stress_max_eigenvalue_over_P_median':float(np.median(ratio)),'ratio_p95':float(np.quantile(ratio,.95)),'ratio_max':float(ratio.max()),'sum_max_eigenvalues_over_sum_thermal_pressure_times_volume':float(ev[mask,1].sum()/((gamma-1)*U[mask].sum())),'anisotropy_energy_weighted':float(np.average(anis,weights=trace[positive])),'cells_ratio_above_one':int(np.sum(ratio>1))}
 rows.append({'subdivisions':N,'regions':regions})
result={'scope':__doc__,'comparison_gamma':gamma,'definition':'Residual stress scale R/V versus P=(gamma-1)U/V; volume cancels in ratio.','rows':rows,'limits':['Residual covariance of diagnostic remapped velocity structure, not molecular thermal stress or an established turbulence model.','No claim that cell-integrated R uniquely determines face traction; a spatial closure is required.','No production pressure or force modified.']}
with (root/'residual_stress_scale.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
