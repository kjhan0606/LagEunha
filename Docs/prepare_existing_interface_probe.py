from pathlib import Path
import json,importlib.util,numpy as np
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_existing_interface_backend_20261018';out.mkdir(exist_ok=False)
spec=importlib.util.spec_from_file_location('rd',base/'run_canonical_av_branches_20261018.py');rd=importlib.util.module_from_spec(spec);spec.loader.exec_module(rd);a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');faces=json.loads((base/'Codex_checkpoint_weight_remap_20261018/certified_face_segments.json').read_text())['faces'];rows=[];gamma=5/3
for face in faces:
 n=np.array(face['normal']);w0,w1=face['normal_speed_coefficients']
 for lo,hi,i,j in face['segments']:
  rhoL,pL=float(a['den'][i]),float(a['pressure'][i]);rhoR,pR=float(a['den'][j]),float(a['pressure'][j]);uL=float(np.array([a['vx'][i],a['vy'][i]])@n);uR=float(np.array([a['vx'][j],a['vy'][j]])@n)
  rows.append([rhoL,pL,uL,np.sqrt(gamma*pL/rhoL),rhoR,pR,uR,np.sqrt(gamma*pR/rhoR),w0+w1*(lo+hi)/2,gamma,0.])
actual=len(rows)
controls=[]
for phase in (0,1):
 for pL,pR,uL,uR in [(1.,1000.,0.,0.),(1e-5,1.,0.,0.),(.01,.01,-1.,1.),(.01,.01,1.,-1.),(1.,1.,0.,0.)]:
  rows.append([1.,pL,uL,np.sqrt(gamma*pL),.01,pR,uR,np.sqrt(gamma*pR/.01),.17,gamma,float(phase)]);controls.append({'index':len(rows)-1,'phase1':phase,'pL':pL,'pR':pR,'uL':uL,'uR':uR})
np.array(rows,dtype=np.float64).tofile(out/'inputs.bin');(out/'inputs.json').write_text(json.dumps({'actual_rows':actual,'total_rows':len(rows),'controls':controls,'source_time':t,'scope':'Cell-constant original donor states at segment midpoint, not fitted-density reconstructed face states.'},indent=2)+'\n')
