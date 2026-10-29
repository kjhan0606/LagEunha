"""Independent geometry and native-log gates for Q2 pressure comparison."""
import hashlib
import json
import re
import socket
import sys
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK,ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
root=Path(sys.argv[1]);summary=Path(sys.argv[2]);out=Path(sys.argv[3])
contract=json.loads((root/'contract.json').read_text());comparison=json.loads(summary.read_text())
original=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/t001_v2/khout.000000.dat')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(original)!='a4d557ae5837e933c8c2c8c335f91e5c97d36f9a3bdb8acd5568b0499140b665':raise ValueError('original IC changed')
a,meta=ordered_global(original);ids=a['u4if']&ID_MASK
r0=np.hypot(a['x']-3,a['y']-3);active=np.flatnonzero(r0<2.)
lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
domain=.5*1024*2.75**2*np.sin(2*np.pi/1024);apothem=2.75*np.cos(np.pi/1024)
gamma=float(np.float32(1.666667));cases=[]
for record in comparison['cases']:
    factor=record['pressure_scale'];d=root/f'x{factor}';path=Path(record['endpoint'])
    if path.parent!=d or sha(path)!=record['endpoint_sha256']:raise ValueError('endpoint changed')
    p,mt=ordered_global(path,ids)
    if abs(mt['t']-contract['target_time'])>1e-9:raise ValueError('time mismatch')
    xy=np.column_stack((p['x']-3.,p['y']-3.));r=np.linalg.norm(xy,axis=1)
    areas,geometry=power_areas(xy,np.zeros(len(a)),active,lookup)
    if geometry['hidden_active'] or np.any(areas<=0):raise ValueError('invalid reconstructed geometry')
    area_max=float(np.max(abs(areas/p['volume'][active].astype(float)-1)))
    vol=p['volume'].astype(float);rho=p['den'].astype(float);m=p['mass'].astype(float)
    pressure=p['pressure'].astype(float);u=p['ie'].astype(float)
    if np.any(vol<=0) or np.any(rho<=0) or np.any(pressure<=0) or np.any(u<=0):raise ValueError('invalid thermodynamics')
    rhoerr=float(np.max(abs(m/(rho*vol)-1)));perr=float(np.max(abs((gamma-1)*u/(vol*pressure)-1)))
    raw=(d/'log').read_text(errors='replace').splitlines()
    endpoint_step=int(path.name.split('.')[1])
    marker=next((k for k,line in enumerate(raw) if line.startswith('Time=') and
                 re.search(rf'\bstep=\s*{endpoint_step}\b',line)),None)
    if marker is None:raise ValueError('missing endpoint step marker')
    log='\n'.join(raw[:marker+1]);lines=[line for line in log.splitlines() if '[RK4E]' in line]
    if not lines:raise ValueError('missing RK4 ledger')
    def value(key):
        match=re.search(rf'\b{key}=([0-9.eE+-]+)',lines[-1]);
        if not match:raise ValueError('missing '+key)
        return float(match.group(1))
    halo=[int(s) for s in re.findall(r'halo_clipped_cells=(\d+)',log)]
    checks={'fixed_polygon':'[FIXED_POLYGON] enabled' in log,
            'fixed_mass_varpress':'[VARPRESS_W0]' in log,
            'inside_polygon':bool(np.all(r<apothem)),
            'area_tiles_polygon':abs(float(np.sum(vol,dtype=np.longdouble))-domain)<1e-5,
            'independent_active_geometry':area_max<1e-5,
            'density_and_eos':max(rhoerr,perr)<1e-6,
            'no_injected_energy':all(value(k)==0 for k in ('floor_cum','de_cum','sfl_cum')),
            'no_geometry_failure':not any(t in log for t in ('[GPU_TESS_INVALID]','[GPU_CLIP_CAPACITY]','[FIXED_POLYGON_CROSSING]')),
            'halo_not_clipped':bool(halo and max(halo)==0),
            'energy_relative_below_1e-3':abs(record['relative_total_energy_change'])<1e-3}
    checks={k:bool(v) for k,v in checks.items()}
    cases.append({'pressure_scale':factor,'endpoint_sha256':record['endpoint_sha256'],
                  'active_sites':len(active),'max_active_area_error':area_max,
                  'max_density_error':rhoerr,'max_eos_error':perr,'max_clipped_halo':max(halo) if halo else None,
                  'checks':checks,'pass':all(checks.values())})
result={'source_comparison_sha256':sha(summary),'cases':cases,'script_sha256':sha(Path(__file__)),
        'pass':all(x['pass'] for x in cases),
        'scope':'Short diagnostic checks only: no claim x4 is original physical IC, no linear stability or t600 proof; all comparisons use own thermal-energy reference.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['pass']:raise SystemExit(1)
