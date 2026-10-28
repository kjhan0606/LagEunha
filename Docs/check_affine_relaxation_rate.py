"""C moment RHS against finite affine pair map and RK4 time refinement."""
import ctypes,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_affine_pair_relaxation_20261018')
from affine_pair_relaxation import relax
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_affine_relaxation_rate_20261018';root.mkdir(exist_ok=False)
header=base/'affine_relaxation_rate.h';(root/header.name).write_bytes(header.read_bytes())
(root/'probe.c').write_text('#include "affine_relaxation_rate.h"\nint probe(const double*a,const double*b,const double*f,const double*t,double g,double*da,double*db){return arr_pair(a,b,f,t,g,da,db);}\n')
subprocess.run(['gcc','-O2','-shared','-fPIC',str(root/'probe.c'),'-o',str(root/'probe.so'),'-lm'],check=True)
lib=ctypes.CDLL(str(root/'probe.so'));ptr=ctypes.POINTER(ctypes.c_double)
lib.probe.argtypes=[ptr,ptr,ptr,ptr,ctypes.c_double,ptr,ptr];lib.probe.restype=ctypes.c_int
def unpack(q):
    m=q[0];c=q[1:3]/m;v=q[6:8]/m
    S=np.array([[q[3],q[4]],[q[4],q[5]]])-m*np.outer(c,c)
    C=q[8:12].reshape(2,2)-m*np.outer(c,v)
    T=np.array([[q[12],q[13]],[q[13],q[14]]])-m*np.outer(v,v)
    G=np.linalg.solve(S,C).T
    return dict(m=m,c=c,v=v,S=S,G=G,R=T-G@S@G.T,U=q[15]-.5*(q[12]+q[14]))
def pack(s,old):
    q=old.copy();m=s['m'];v=s['v'];G=s['G'];S=s['S']
    q[6:8]=m*v;q[8:12]=(S@G.T+m*np.outer(s['c'],v)).ravel()
    Pi=G@S@G.T+s['R']+m*np.outer(v,v)
    q[12:15]=[Pi[0,0],Pi[0,1],Pi[1,1]];q[15]=s['U']+.5*np.trace(Pi)
    return q
path=base/'Codex_checkpoint_weight_remap_20261018/state_w1_q8.npz'
with np.load(path) as data:q=data['state']
a=unpack(q[0]);c=q[:,1:3]/q[:,0,None];dist=np.linalg.norm(c-a['c'],axis=1);dist[0]=np.inf
j=int(np.argmin(dist));b=unpack(q[j]);delta=b['c']-a['c'];t=np.array([-delta[1],delta[0]])/np.linalg.norm(delta)
f=np.ascontiguousarray(.5*(a['c']+b['c'])+.2*np.linalg.norm(delta)*t);t=np.ascontiguousarray(t)
def rhs(y):
    out=np.empty_like(y);code=lib.probe(y[0].ctypes.data_as(ptr),y[1].ctypes.data_as(ptr),f.ctypes.data_as(ptr),t.ctypes.data_as(ptr),.5,out[0].ctypes.data_as(ptr),out[1].ctypes.data_as(ptr))
    assert code==0
    return out
y0=np.ascontiguousarray([q[0],q[j]]);r=rhs(y0)
_,_,info=relax(a,b,f,t,.5,1.);duration=.4/(.5*info['mobility'])
ea,eb,_=relax(a,b,f,t,.5,duration);exact=np.array([pack(ea,y0[0]),pack(eb,y0[1])])
rows=[]
for steps in [4,8,16]:
    y=y0.copy();dt=duration/steps
    for k in range(steps):
        k1=rhs(y);k2=rhs(y+.5*dt*k1);k3=rhs(y+.5*dt*k2);k4=rhs(y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
    rows.append(dict(steps=steps,dt=dt,max_error_vs_exact_map=float(abs(y-exact).max()),
                     total_energy_change=float(sum(y[:,15])-sum(y0[:,15])),
                     total_angular_change=float(sum(y[:,9]-y[:,10])-sum(y0[:,9]-y0[:,10]))))
errors=[x['max_error_vs_exact_map'] for x in rows];assert errors[2]<errors[1]<errors[0]
assert max(abs(x['total_energy_change']) for x in rows)<1e-15
out=dict(scope=__doc__,duration=duration,pair_indices=[0,j],rows=rows,
         error_ratios=[errors[k]/errors[k+1] for k in range(2)],
         instantaneous_energy_rate_sum=float(sum(r[:,15])),instantaneous_angular_rate_sum=float(sum(r[:,9]-r[:,10])),
         thermal_rate_per_cell=(r[:,15]-.5*(r[:,12]+r[:,14])).tolist(),
         header_sha256=hashlib.sha256(header.read_bytes()).hexdigest(),
         limits=['One frozen pair, host C only; no complete mesh or moving Euler RK4 claim.',
                 'RK4 positive heat/state is not guaranteed for arbitrary dt.',
                 'Not connected to production Kepler or MPI.'])
(root/'analysis.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
