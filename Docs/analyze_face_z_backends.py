"""Analyze the isolated FACE_Z probe on the compute host."""
import json
from pathlib import Path
import numpy as np

p = Path('/gpfs/kjhan/LagForce/Codex_active_face_z_backend_20261018')
a = np.fromfile(p/'inputs.bin').reshape(-1, 13)
c = np.fromfile(p/'cpu.bin').reshape(-1, 6)
g = np.fromfile(p/'gpu.bin').reshape(-1, 6)
r = np.load(p/'recorded_ZQ.npy')
mm = np.where(a[:,9]*a[:,10] > 0,
              np.sign(a[:,9])*np.minimum(abs(a[:,9]),abs(a[:,10])), 0)
fi, fj = np.maximum(a[:,7],1e-12), np.maximum(a[:,8],1e-12)
rl, rr = a[:,5]*fi, a[:,6]*fj
cm = np.maximum(np.sqrt(a[:,12]*a[:,3]*fi**a[:,12]/rl),
                np.sqrt(a[:,12]*a[:,4]*fj**a[:,12]/rr))
mu = rl*rr/(rl+rr)
ah = (a[:,12]+1)*abs(mm)/8
z = mu*np.maximum(cm+np.maximum(mm,0), ah+np.sqrt(ah**2+cm**2))
q = z*mm
ref = np.column_stack((z,q,-a[:,0]*q*a[:,1],-a[:,0]*q*a[:,2],
                       -a[:,0]*q*a[:,11],.5*a[:,0]*q*a[:,10]))
result = {
    'rows':len(a), 'output_columns':['Zq','Q','FQx','FQy','EQ','heat_i'],
    'cpu_gpu_max_abs_by_column':np.max(abs(c-g),axis=0).tolist(),
    'numpy_reference_max_abs_by_column':np.max(abs(c-ref),axis=0).tolist(),
    'recorded_ZQ_relative_l2':(np.linalg.norm(c[:,:2]-r,axis=0)/np.linalg.norm(r,axis=0)).tolist(),
    'negative_directed_heat_count':int(np.count_nonzero(c[:,5]<0)),
    'scope':'Recorded rounded upstream inputs; isolated dissipative algebra only. No full force, MPI, RK4, or trajectory validation.',
}
with (p/'analysis.json').open('x') as f:
    json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
