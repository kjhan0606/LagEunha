"""C spin ledger fed by actual affine ALE patch fluxes and volume sources.

Independent differences rebuild moving power cells and integrate a stationary
Euler shear field. Small smooth patch, not Kepler or shock validation.
"""
import ctypes,json,hashlib
from pathlib import Path
import numpy as np
from affine_moment_patch import affine_case
root=Path('/gpfs/kjhan/LagForce/Codex_moment_spin_ledger_20261018')
lib=ctypes.CDLL(str(root/'probe.so'));ptr=ctypes.POINTER(ctypes.c_double)
lib.probe.argtypes=[ptr,ptr,ptr];lib.probe.restype=ctypes.c_int
def pack(raw):
    assert np.max(abs(raw[:,4]-raw[:,5]))<1e-12
    assert np.max(abs(raw[:,14]-raw[:,15]))<1e-12
    q=np.zeros((len(raw),17));q[:,:16]=raw[:,[0,1,2,3,4,6,7,8,9,10,11,12,13,14,16,17]]
    return q
def calculate(q,d):
    result=[]
    for a,b in zip(q,d):
        z=np.empty(8)
        assert lib.probe(a.ctypes.data_as(ptr),b.ctypes.data_as(ptr),z.ctypes.data_as(ptr))==0
        result.append(z)
    return np.array(result)
cases=[]
for n in [4,8]:
    captured={};baseline=affine_case(n,capture=captured)
    q=pack(captured['raw']);d=pack(captured['rate']);o=calculate(q,d)
    errors=[]
    for eps,plus,minus in captured['finite_states']:
        a=pack(plus);b=pack(minus)
        za=calculate(a,np.zeros_like(a));zb=calculate(b,np.zeros_like(b))
        derivative=(za-zb)/(2*eps)
        # omega_dot has its own chain rule, not an inertia-rate column.
        expected=np.column_stack((o[:,4],o[:,5],(o[:,5]-o[:,2]*o[:,4])/o[:,0],o[:,6]))
        errors.append(dict(eps=eps,max_abs_errors_I_L_omega_Krot=np.max(abs(derivative[:,:4]-expected),axis=0).tolist()))
    assert errors[1]['max_abs_errors_I_L_omega_Krot'][3]<errors[0]['max_abs_errors_I_L_omega_Krot'][3]
    assert baseline['pressure_rate_max']<1e-10
    cases.append(dict(n=n,cells=n*n,baseline=baseline,spin_finite_difference=errors,
        max_Krot_rate=float(abs(o[:,6]).max()),sum_Krot_rate=float(o[:,6].sum()),
        max_thermal_rate=float(abs(o[:,7]).max()),
        max_cell_mass_rate=float(abs(d[:,0]).max()),
        max_internal_angular_rate=float(abs(o[:,5]).max())))
print(json.dumps(dict(scope=__doc__,cases=cases,
    header_sha256=hashlib.sha256((root/'moment_spin_ledger.h').read_bytes()).hexdigest(),
    limits=['Fixed net cell mass via weight-rate solve, but local face mass flux is nonzero.',
    'Thermal and spin rates must not be conflated. No relaxation law added.',
    'Only host C ledger tested; full GPU/MPI evolution remains incomplete.']),indent=2))
