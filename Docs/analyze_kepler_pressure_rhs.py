"""Material-bin pressure-rate decomposition from a pre-HG canonical RHS.

Uses geometric volume rate in the instantaneous chain rule. Does not replace
a finite-step GCL check or include split HG heating.
"""
import argparse
import hashlib
import json
import re
import struct
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--run',required=True,type=Path)
p.add_argument('--initial',required=True,type=Path)
p.add_argument('--layout',required=True,type=Path)
p.add_argument('--output',required=True,type=Path)
a=p.parse_args()
layout=json.loads(a.layout.read_text())
raw=a.initial.read_bytes()
magic,n,size,stress=struct.unpack_from('<iiii',raw)
assert magic==-1 and stress==1 and size==layout['size']
assert len(raw)==16+n*size+52
assert struct.unpack_from('<d',raw,16+n*size)[0]==0
names=['u4if','x','y','mass','pressure','den']
dtype=np.dtype({'names':names,'formats':['<u8']+['<f8']*5,
                'offsets':[layout[k] for k in names],'itemsize':size})
initial=np.frombuffer(raw,dtype=dtype,count=n,offset=16)
ids=initial['u4if'] & layout['id_mask']
assert np.array_equal(np.sort(ids),np.arange(n))
initial=initial[np.argsort(ids)]
c=np.loadtxt(a.run/'rhs.dat');c=c[np.argsort(c[:,0])]
assert c.shape==(n,25) and np.array_equal(c[:,0],np.arange(n))
assert np.isfinite(c).all() and np.all(c[:,6]>0)
assert np.allclose(c[:,7],initial['mass'],rtol=1e-12,atol=0)
gamma=float(np.float32(float(re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)',(a.run/'params.dat').read_text())[1])))
p0=float(initial['pressure'][0]);assert np.all(initial['pressure']==p0)
r0=np.hypot(initial['x']-3,initial['y']-3)
r=np.hypot(c[:,1]-3,c[:,2]-3)
du=c[:,10]-np.sum(c[:,3:5]*c[:,8:10],axis=1)
compression=-gamma*c[:,5]*c[:,12]/c[:,6]/p0
heating=(gamma-1)*c[:,11]/c[:,6]/p0
direct=((gamma-1)*du-c[:,5]*c[:,12])/c[:,6]/p0
residual=direct-compression-heating
regions={'inner_background':(initial['den']<.5)&(r0<=.5),
         'disk_inner_edge':(initial['den']>.5)&(r0<.6)}
for lo,hi in [(0,.125),(.125,.25),(.25,.375),(.375,.5)]:
    regions[f'inner_r0_{lo:g}_{hi:g}']=(initial['den']<.5)&(r0>=lo)&(r0<hi)
out={}
for name,sel in regions.items():
    assert np.all(r[sel]<2.5)
    w=initial['mass'][sel]
    mean=lambda x: float(np.average(x[sel],weights=w))
    fluct=c[:,5]/p0-mean(c[:,5]/p0)
    error=float(np.abs(residual[sel]).sum()/np.maximum(np.abs(direct[sel]).sum(),1e-300))
    assert error<1e-10,(name,error)
    out[name]={'n':int(sel.sum()),'mean_pressure_rate_over_p0':mean(direct),
               'compression_mean_rate':mean(compression),'TCF_mean_rate':mean(heating),
               'variance_growth_rate':2*mean(fluct*direct),
               'compression_variance_growth':2*mean(fluct*compression),
               'TCF_variance_growth':2*mean(fluct*heating),
               'chain_rule_relative_L1_residual':error}
result={'source':str(a.run),'rhs_sha256':hashlib.sha256((a.run/'rhs.dat').read_bytes()).hexdigest(),
        'initial_sha256':hashlib.sha256(raw).hexdigest(),'gamma':gamma,'pressure0':p0,
        'equation':'Pdot/P0 = -gamma*(P/P0)*Vdot/V + (gamma-1)*TCF_heat/(P0*V)',
        'scope':'Instantaneous pre-HG RHS, fixed initial material membership; split HG excluded; geometric Vdot chain rule, not a finite-step pressure budget.',
        'regions':out}
with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
