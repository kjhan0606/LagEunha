"""Generator velocity and internal kinetic split after saved hybrid remap paths."""
from pathlib import Path
import ast,json,importlib.util
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_hybrid_affine_path_20261018'
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
rd=load('rd',base/'run_canonical_av_branches_20261018.py');a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));vel=np.column_stack((a['vx'],a['vy']))
tree=ast.parse((base/'audit_remap_true_spatial.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='recover'],type_ignores=[]),'recover','exec'),globals())
rows=[]
for N in (4,8):
 saved=np.load(root/f'steps{N}.npz');q=saved['state'];ids=saved['ids'];m=q[:,0];c=q[:,1:3]/m[:,None];v=q[:,6:8]/m[:,None];S,C,G,R=recover(q);vg=v+np.einsum('nij,nj->ni',G,x[ids]-c);r=np.linalg.norm(x[ids],axis=1);er=x[ids]/np.maximum(r,1e-300)[:,None];et=np.column_stack((-er[:,1],er[:,0]));aff=.5*np.trace(G@S@G.transpose(0,2,1),axis1=1,axis2=2);res=.5*np.trace(R,axis1=1,axis2=2);U=q[:,15]-.5*(q[:,12]+q[:,14]);regions={}
 for name,mask in [('patch',r<.7),('inner',r<.5),('edge',(r>=.5)&(r<.6)),('outer',(r>=.6)&(r<.7))]:
  weight=a['mass'][ids][mask];dmean=v-vel[ids];dgen=vg-vel[ids]
  def rms(z):return float(np.sqrt(np.average(z[mask]**2,weights=weight)))
  regions[name]={'cells':int(mask.sum()),'mean_velocity_change_RMS':rms(np.linalg.norm(dmean,axis=1)),'generator_velocity_change_RMS':rms(np.linalg.norm(dgen,axis=1)),'generator_radial_change_RMS':rms(np.sum(dgen*er,axis=1)),'generator_tangential_change_RMS':rms(np.sum(dgen*et,axis=1)),'affine_K_over_original_U':float(aff[mask].sum()/a['ie'][ids][mask].sum()),'residual_K_over_original_U':float(res[mask].sum()/a['ie'][ids][mask].sum()),'thermal_change_over_original_U':float((U[mask].sum()-a['ie'][ids][mask].sum())/a['ie'][ids][mask].sum())}
 rows.append({'subdivisions':N,'regions':regions})
result={'scope':__doc__,'source_time':t,'rows':rows,'limits':['Velocity changes relative to the original t2 point values, not errors against a known evolved solution.','Generator velocity is diagnostic recovery, not adopted generator motion or steering.','Original uniform lift has constant velocity; no added initial affine energy.','Residual covariance is nonthermal. No force or physical time step.']}
with (root/'generator_recovery.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
