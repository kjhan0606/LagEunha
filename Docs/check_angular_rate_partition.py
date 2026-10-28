"""Exercise the actual diagnostic function, provided through stdin AST."""
import ast,json,sys
import numpy as np
source=sys.stdin.read()
node=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='angular_rate_partition')
exec(compile(ast.Module(body=[node],type_ignores=[]),'<actual diagnostic>','exec'))
m=np.array([2.,3.]);raw=np.zeros((2,18));rates=np.zeros_like(raw)
raw[:,1:3]=[[4.,6.],[-3.,6.]];raw[:,7:9]=[[1.,-2.],[4.,5.]]
rates[:,0]=[.02,0.];rates[:,1:3]=[[.4,-.1],[0.,0.]]
rates[:,7:9]=[[.8,.3],[0.,0.]];rates[:,10]=[.7,0.];rates[:,11]=[-.2,0.]
H=np.array([[[0.,.4],[-.1,0.]],[[0.,8.],[1.,0.]]])
r=angular_rate_partition(raw,rates,m,H,np.array([False,True]))
def ell(a,m):return a[:,10]-a[:,11]-(a[:,1]*a[:,8]-a[:,2]*a[:,7])/m
h=1e-4
fd=(ell(raw+h*rates,m+h*rates[:,0])-ell(raw-h*rates,m-h*rates[:,0]))/(2*h)
error=float(np.max(abs(fd-r['internal_angular_rate'])))
assert error<1e-9
assert r['pressure_internal_torque'][0]==.5
assert all(v[1]==0 for v in r.values())
assert np.allclose(r['internal_angular_rate'],r['pressure_internal_torque']+r['transport_internal_angular_rate'],rtol=0,atol=1e-15)
print(json.dumps(dict(passed=True,finite_difference_error=error,values={k:v.tolist() for k,v in r.items()},limits='Algebraic diagnostic fixture only, no actual RHS evaluation.'),indent=2))
