from pathlib import Path
import json,struct,hashlib,numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_t2_volume_derivative_20261017');source=Path('/gpfs/kjhan/LagForce/Codex_t2_onset_rhs_20261017/circular/khout.000184.dat');original=source.read_bytes();layout=json.loads(Path('/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json').read_text());magic,n,size,stress=struct.unpack_from('<iiii',original);assert magic==-1 and size==layout['size']
names=['x','y','vx','vy'];dtype=np.dtype({'names':names,'formats':['<f8']*4,'offsets':[layout[k] for k in names],'itemsize':size});a=np.frombuffer(original,dtype=dtype,count=n,offset=16);x=a['x']-3;y=a['y']-3;r2=x*x+y*y;omega=(x*a['vy']-y*a['vx'])/r2
reports=[]
for e in json.loads((root/'manifest.json').read_text()):
 raw=bytearray(original);b=np.ndarray((n,),dtype=dtype,buffer=raw,offset=16);theta=e['signed_dt']*omega;c=np.cos(theta);s=np.sin(theta);xx=c*x-s*y;yy=s*x+c*y;b['x']=np.mod(xx+3,6);b['y']=np.mod(yy+3,6);b['vx']=-omega*yy;b['vy']=omega*xx
 changed=np.flatnonzero(np.frombuffer(original,dtype=np.uint8)!=np.frombuffer(raw,dtype=np.uint8));off=(changed-16)%size;allowed=np.zeros(len(changed),bool)
 for field in names:allowed|=(off>=layout[field])&(off<layout[field]+8)
 assert np.all(allowed&(changed>=16)&(changed<16+n*size))
 p=root/e['case']/'khout.000184.dat'
 with p.open('xb') as f:f.write(raw)
 reports.append(dict(e,sha256=hashlib.sha256(raw).hexdigest(),only_position_velocity_bytes_changed=True))
(root/'inputs.json').write_text(json.dumps({'source':str(source),'source_sha256':hashlib.sha256(original).hexdigest(),'scope':'Exact circular rotation of diagnostic copies; all other particle fields and time header retained. Not physical evolution. Analysis restricted away from periodic/frozen boundary.','cases':reports},indent=2))
