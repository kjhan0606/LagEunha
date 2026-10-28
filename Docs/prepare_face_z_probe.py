from pathlib import Path
import numpy as np,json
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_active_face_z_backend_20261018';out.mkdir(exist_ok=False);root=base/'Codex_TCF_rhs_identity_20261017/t6';f=np.loadtxt(root/'rhs.dat.faces');c=np.loadtxt(root/'rhs.dat');c=c[np.argsort(c[:,0])];i=f[:,0].astype(int);j=f[:,1].astype(int);mean=.5*np.sum((c[i,3:5]+c[j,3:5])*f[:,6:8],axis=1)
a=np.column_stack((f[:,4],f[:,6:8],c[i,5],c[j,5],c[i,7]/c[i,6],c[j,7]/c[j,6],f[:,16:19],f[:,9],mean,np.full(len(f),1.6666669845581055)));assert a.shape[1]==13;a.tofile(out/'inputs.bin');np.save(out/'recorded_ZQ.npy',f[:,[15,14]])
(out/'input_scope.json').write_text(json.dumps({'rows':len(f),'source':str(root),'gamma':1.6666669845581055,'scope':'Recorded rounded face qin,fi,fj,du plus saved owner rho/P and mean normal velocity; upstream reconstruction not reproduced.'},indent=2)+'\n')
