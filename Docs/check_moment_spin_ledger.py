"""C chain-rule check on existing positive remap states and a specified path.

The path is the chord between two saved remaps; it is not a physical RHS.
"""
import ctypes,hashlib,json,subprocess
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
root=base/'Codex_moment_spin_ledger_20261018'
root.mkdir(exist_ok=False)
header=base/'moment_spin_ledger.h'
(root/'moment_spin_ledger.h').write_bytes(header.read_bytes())
(root/'probe.c').write_text('#include "moment_spin_ledger.h"\nint probe(const double*q,const double*d,double*out){return msl_rate(q,d,out);}\n')
cmd=['gcc','-O2','-shared','-fPIC',str(root/'probe.c'),'-o',str(root/'probe.so'),'-lm']
subprocess.run(cmd,check=True)
lib=ctypes.CDLL(str(root/'probe.so'));fn=lib.probe
ptr=ctypes.POINTER(ctypes.c_double);fn.argtypes=[ptr,ptr,ptr];fn.restype=ctypes.c_int
def call(q,d):
    o=np.full(8,12345.,dtype=float)
    code=fn(q.ctypes.data_as(ptr),d.ctypes.data_as(ptr),o.ctypes.data_as(ptr))
    return code,o
paths=[base/'Codex_checkpoint_weight_remap_20261018'/f'state_w{f}_q8.npz' for f in ['0.25','1']]
states=[];ids=[]
for p in paths:
    with np.load(p) as a:states.append(a['state']);ids.append(a['ids'])
assert np.array_equal(*ids)
q=np.ascontiguousarray((states[0]+states[1])*.5)
d=np.ascontiguousarray(states[1]-states[0])
def reference(q):
    # Long double raw-moment evaluation for an independent difference check.
    a=q.astype(np.longdouble);m=a[:,0];cx=a[:,1]/m;cy=a[:,2]/m
    I=a[:,3]+a[:,5]-cx*a[:,1]-cy*a[:,2]
    L=a[:,9]-a[:,10]-cx*a[:,7]+cy*a[:,6]
    assert np.all(I>0)
    return .5*L*L/I
outputs=[]
for a,b in zip(q,d):
    code,o=call(a,b);assert code==0;outputs.append(o)
o=np.array(outputs);errors=[]
for h in [.002,.001,.0005]:
    finite=(reference(q.astype(np.longdouble)+h*d)-reference(q.astype(np.longdouble)-h*d))/(2*h)
    errors.append(float(np.max(abs(finite-o[:,6]))))
assert errors[-1]<errors[0]
# Mass-changing paths must include mdot terms; compare deliberately omitted ones.
no_mass=d.copy();no_mass[:,0]=0
omitted=np.array([call(a,b)[1] for a,b in zip(q,no_mass)])
guards=[]
for mode in ['zero_mass','negative_inertia','nan']:
    a=q[0].copy()
    if mode=='zero_mass':a[0]=0
    elif mode=='negative_inertia':a[3]=a[5]=-1
    else:a[7]=np.nan
    code,z=call(a,d[0]);assert code!=0 and np.all(z==12345.)
    guards.append(dict(case=mode,code=code,output_unchanged=True))
out=dict(scope=__doc__,cells=len(q),header_sha256=hashlib.sha256(header.read_bytes()).hexdigest(),
    input_sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in paths],compiler_command=cmd,
    finite_difference_max_errors=errors,successive_error_ratios=[errors[i]/errors[i+1] for i in range(2)],
    max_absolute_mass_rate=float(abs(d[:,0]).max()),
    max_change_Krot_rate_if_mass_terms_omitted=float(abs(o[:,6]-omitted[:,6]).max()),guards=guards,
    limits=['C host diagnostic only; CUDA execution and actual ALE RHS not checked here.',
            'No production energy, force, weight law or relaxation update was added.'])
(root/'analysis.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
