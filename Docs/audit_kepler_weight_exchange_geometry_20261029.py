"""Check weight-only mass flux against independent power-cell area differences."""
import ctypes
import hashlib
import json
import socket
import sys
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
prior=json.loads((docs/'KEPLER_Q2_NORMAL_TRACTION_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
path=root/'khout.011058.dat';initial=Path(late['initial'])
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
a,ma=ordered_global(initial);p,mt=ordered_global(path,a['u4if']&ID_MASK)
ref=next(v for v in prior['cases'] if v['time']==mt['t'])
if sha(initial)!=late['initial_sha256'] or sha(path)!=ref['source_sha256'] or ma['t']!=0 or np.any(p['w2']!=0):raise ValueError('changed input')
xy0=np.column_stack((a['x']-3,a['y']-3));r0=np.linalg.norm(xy0,axis=1)
selected=np.flatnonzero((r0>=.5)&(r0<.75));flag=np.zeros(len(a),bool);flag[selected]=True
lookup=np.full(len(a),-1,int);lookup[selected]=np.arange(len(selected))
xy=np.column_stack((p['x']-3,p['y']-3));zero=np.zeros(len(a))
area,info,tri,centers=power_areas(xy,zero,selected,lookup,return_geometry=True)
if info['hidden_active'] or np.max(abs(area/p['volume'][selected]-1))>1e-5:raise ValueError('bad original geometry')
i,j,c,length=finite_faces_touching_active(tri,centers,flag)
dist=np.linalg.norm(xy[j]-xy[i],axis=1)
rates=np.random.default_rng(20261029).normal(0,1e-4,len(a))
lib=ctypes.CDLL(str(Path(sys.argv[1]).resolve()));fn=lib.weight_exchange_batch
arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
fn.argtypes=[ctypes.c_int,arr,arr,arr,arr];fn.restype=ctypes.c_int
di=np.ascontiguousarray(rates[i]);dj=np.ascontiguousarray(rates[j]);dist=np.ascontiguousarray(dist);flux=np.empty(len(i))
if fn(len(i),di,dj,dist,flux)!=1:raise ValueError('candidate flux rejected')
rate=np.zeros(len(selected));ii=flag[i];jj=flag[j]
np.add.at(rate,lookup[i[ii]],-length[ii]*flux[ii]);np.add.at(rate,lookup[j[jj]],length[jj]*flux[jj])
rows=[]
for h in (1e-4,5e-5,2.5e-5):
    plus,meta_plus=power_areas(xy,h*rates,selected,lookup)
    minus,meta_minus=power_areas(xy,-h*rates,selected,lookup)
    if meta_plus['hidden_active'] or meta_minus['hidden_active']:raise ValueError('perturbed cell hidden')
    diff=(plus-minus)/(2*h);error=diff-rate
    rows.append({'step':h,'rms_area_derivative':float(np.sqrt(np.mean(diff**2))),
                 'rms_difference':float(np.sqrt(np.mean(error**2))),
                 'max_abs_difference':float(np.max(abs(error))),
                 'relative_rms_difference':float(np.linalg.norm(error)/np.linalg.norm(diff))})
result={'snapshot_time':mt['t'],'snapshot_sha256':ref['source_sha256'],'initial_sha256':late['initial_sha256'],
        'sites':len(selected),'faces':len(i),'weight_rate_seed':20261029,'weight_rate_std':1e-4,
        'steps':rows,'bridge_sha256':sha(Path(sys.argv[1])),'script_sha256':sha(Path(__file__)),
        'scope':'Frozen Q2 t47.5 bounded W0 geometry. Unit-density weight-only flux from actual C helper compared to centered differences of independently rebuilt power-cell areas at fixed generators. This validates infinitesimal geometry/sign only; no fluid evolution, face-state selection, full ALE flux, finite-time transport, topology-change stability, energy work or MPI.'}
out=docs/'KEPLER_WEIGHT_EXCHANGE_GEOMETRY_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
