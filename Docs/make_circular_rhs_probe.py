"""Create a diagnostic-only snapshot copy, replacing vx/vy and nothing else.

Never use this artificial state as a production restart. It is intended for
one force evaluation terminated by GFS_TCF_DUMP_CALL=GFS_TCF_DUMP_LAST=1.
"""
import argparse
import hashlib
import json
import struct
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser()
ap.add_argument('source')
ap.add_argument('destination')
ap.add_argument('--layout',required=True)
args=ap.parse_args()
layout=json.loads(Path(args.layout).read_text())
raw=Path(args.source).read_bytes()
magic,n,size,stress=struct.unpack_from('<iiii',raw)
assert magic==-1 and stress==1 and size==layout['size']
assert len(raw)==16+n*size+52
buf=bytearray(raw)
dt=np.dtype({'names':['x','y','vx','vy'], 'formats':['<f8']*4,
             'offsets':[layout[k] for k in ['x','y','vx','vy']], 'itemsize':size})
a=np.frombuffer(buf,dtype=dt,count=n,offset=16)
x,y=a['x']-3,a['y']-3
r=np.sqrt(x*x+y*y); u=r/.5
mass=np.ones_like(r)
lo=u<=.5; mid=(u>.5)&(u<1)
q=u[lo]; mass[lo]=q**3*(32/3-192/5*q*q+32*q**3)
q=u[mid]; mass[mid]=-1/15+q**3*(64/3-48*q+192/5*q*q-32/3*q**3)
omega=np.sqrt(np.maximum(mass,0)/np.maximum(r,1e-30)**3)
a['vx']=-y*omega; a['vy']=x*omega
assert np.isfinite(a['vx']).all() and np.isfinite(a['vy']).all()
# Restore just the permitted bytes and require equality to the whole input.
check=bytearray(buf)
for i in range(n):
    for k in ['vx','vy']:
        j=16+i*size+layout[k]
        check[j:j+8]=raw[j:j+8]
assert check==raw
dest=Path(args.destination)
with dest.open('xb') as f: f.write(buf)
manifest={'source':args.source,'destination':str(dest),'particles':n,
          'source_sha256':hashlib.sha256(raw).hexdigest(),
          'probe_sha256':hashlib.sha256(buf).hexdigest(),
          'only_vx_vy_changed':True,'center':[3,3],'softening':'spline','rc':.5,
          'scope':'Artificial velocity field for a single RHS diagnostic; not an evolved solution.'}
dest.with_suffix('.probe.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest))
