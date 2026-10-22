"""Compare continuous-field initialization integrals with legacy point state.
This is a read-only audit, not a conservative conversion or restart method.
"""
from pathlib import Path
import argparse,json,struct,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from audit_snapshot import read_snapshot
p=argparse.ArgumentParser();p.add_argument('--version',default='v4');args=p.parse_args()
dest=ROOT/f'rk4_diagnostics/enriched_{args.version}'
initial=ROOT/'core_refinement/fine_integrated/khout.000001.dat'
a,meta=read_snapshot(initial);assert meta['t']==0
path=dest/'fine_fresh/enriched_project_r0.bin'
with path.open('rb') as f:h=struct.unpack('<8Q',f.read(64))
dtype=np.dtype(dict(names=['id','area','integral'],formats=['<u8','<f8',('<f8',4)],offsets=[0,24,1360],itemsize=h[5]))
b=np.fromfile(path,dtype=dtype,offset=64,count=h[2]);b=b[np.argsort(b['id'])];a=a[np.argsort(a['u4if'])]
assert np.array_equal(a['u4if'],b['id'])
old=np.column_stack((a['mass'],a['mass']*a['vx'],a['mass']*a['vy'],a['ie']+.5*a['mass']*(a['vx']**2+a['vy']**2)))
delta=b['integral']-old
def sums(v):return [float(x) for x in np.sum(v,axis=0,dtype=np.longdouble)]
thermal=float(np.sum(a['ie'],dtype=np.longdouble))
gamma=float(np.float32(1.666667))
continuous_thermal=float(np.sum(b['area'],dtype=np.longdouble))*1.e-6/(gamma-1)
out=dict(scope='fresh continuous-field integrals versus legacy point-lumped cell state; NOT a conservative conversion or restart',
    initial=str(initial),legacy_totals=sums(old),continuous_totals=sums(b['integral']),total_delta=sums(delta),
    legacy_thermal_total=thermal,total_E_delta_over_legacy_thermal=float(np.sum(delta[:,3],dtype=np.longdouble))/thermal,
    gamma=gamma,continuous_thermal_total=continuous_thermal,thermal_delta=continuous_thermal-thermal,
    kinetic_delta=float(np.sum(delta[:,3],dtype=np.longdouble))-(continuous_thermal-thermal),
    max_cell_relative_mass_delta=float(np.max(abs(b['integral'][:,0]/old[:,0]-1))),
    classification='Projection conserves its continuous-field integrals. It does not preserve legacy cell momentum/E; disclose and settle fresh-initialization contract before evolution. No input or checkpoint has been changed.')
(dest/'initial_ledger.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
