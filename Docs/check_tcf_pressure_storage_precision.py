"""Isolate GPU-style pressure storage rounding on preserved actual face inputs.

Both evaluations use float64 arithmetic and identical geometry, velocity, Q.
Only Pi/Pj are rounded through float32. This is not the unsupported GPU TCF
path, not a QBP replay, and not a complete mixed-precision error estimate.
"""
import json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_tcf_full_face_probe_20261018')
a=np.fromfile(root/'inputs.bin',dtype=np.float64).reshape(-1,13)
owners=np.fromfile(root/'owners.bin',dtype=np.int64)
cpu=np.fromfile(root/'cpu.bin',dtype=np.float64).reshape(-1,6)
assert len(a)==len(owners)==len(cpu)
def evaluate(a):
    A,d,ex,ey,tx,ty,Pi,Pj,vix,viy,vjx,vjy,Q=a.T
    ui=vix*ex+viy*ey;uj=vjx*ex+vjy*ey;Pb=.5*(Pi+Pj)
    out=np.zeros((len(a),6))
    out[:,0]=-A*(Pb+Q)*ex-A*(Pj-Pi)*tx/d
    out[:,1]=-A*(Pb+Q)*ey-A*(Pj-Pi)*ty/d
    out[:,2]=-A*(.5*(Pj*ui+Pi*uj)+(Pj*(vix*tx+viy*ty)-Pi*(vjx*tx+vjy*ty))/d+.5*Q*(ui+uj))
    out[:,3]=A*(.5*(ui+uj)+((vix-vjx)*tx+(viy-vjy)*ty)/d)
    out[:,4]=.5*A*Q*(ui-uj)
    return out
reference=evaluate(a);rounded=a.copy();rounded[:,6:8]=a[:,6:8].astype(np.float32).astype(np.float64)
changed=evaluate(rounded)
assert np.max(abs(reference-cpu))<1e-12
n=int(owners.max())+1;totals=np.zeros((n,6));rounded_totals=np.zeros_like(totals)
np.add.at(totals,owners,reference);np.add.at(rounded_totals,owners,changed)
selected=np.unique(owners)
def metrics(ref,value):
    return {'relative_L2':float(np.linalg.norm(value-ref)/max(np.linalg.norm(ref),1e-300)),
            'max_absolute':float(np.max(abs(value-ref)))}
out={'scope':__doc__,'faces':len(a),'owners':len(selected),
    'pressure_storage':metrics(a[:,6:8],rounded[:,6:8]),
    'independent_FP64_vs_saved_CPU_max_abs':float(np.max(abs(reference-cpu))),
    'face_force':metrics(reference[:,:2],changed[:,:2]),
    'summed_cell_force':metrics(totals[selected,:2],rounded_totals[selected,:2]),
    'face_energy_work':metrics(reference[:,2],changed[:,2]),
    'summed_cell_energy_work':metrics(totals[selected,2],rounded_totals[selected,2]),
    'volume_and_heat_unchanged':bool(np.array_equal(reference[:,3:],changed[:,3:])),
    'limitations':'Pressure quantization only. Frozen rounded Q means density/sound-speed/QBP sensitivity is excluded. Full GPU equations remain unimplemented. An observed rounding difference alone does not demonstrate orbital instability.'}
print(json.dumps(out,indent=2))
