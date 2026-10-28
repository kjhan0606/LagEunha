"""Explicit CPU batch API; installed SciPy's public Cython Brent ABI."""
import ctypes
import os
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize.cython_optimize import _zeros

_capsule=_zeros.__pyx_capi__['brentq']
_name=ctypes.pythonapi.PyCapsule_GetName
_name.argtypes=[ctypes.py_object];_name.restype=ctypes.c_char_p
_signature=_name(_capsule)
if b'callback_type' not in _signature or b'zeros_full_output' not in _signature:
    raise RuntimeError('Unsupported SciPy Cython Brent ABI')
_get=ctypes.pythonapi.PyCapsule_GetPointer
_get.argtypes=[ctypes.py_object,ctypes.c_char_p];_get.restype=ctypes.c_void_p
_brent=_get(_capsule,_signature)
_path=Path(os.environ['LAGFORCE_NATIVE_EULER_STAR_LIB']).resolve(strict=True)
_lib=ctypes.CDLL(str(_path))
_array=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
_lib.euler_star_batch.argtypes=[ctypes.c_size_t,_array,_array,ctypes.c_void_p]
_lib.euler_star_batch.restype=None

def provenance():
    return dict(library=str(_path),scipy_version=scipy.__version__,brent_ABI=_signature.decode(),
                xtol=float(np.finfo(float).tiny),rtol=float(8*np.finfo(float).eps),maxiter=200)

def star_batch(inputs,return_failures=False):
    """Rows [rhoL,uL,pL,rhoR,uR,pR,gamma]; output p/speed + root diagnostics.

    Status0 success,1 invalid,2 vacuum,3 bracket overflow,4 Brent failure,
    5 unresolved positive pressure,6 nonfinite normal-range calculation.
    The error envelopes assume normal arithmetic and usual libm accuracy;
    they are not interval-certified bounds for extreme under/overflow.
    """
    inputs=np.ascontiguousarray(inputs,dtype=float)
    if inputs.ndim!=2 or inputs.shape[1]!=7:raise ValueError('expected (records,7) states')
    result=np.empty((len(inputs),13),dtype=float)
    _lib.euler_star_batch(len(inputs),inputs,result,_brent)
    failed=np.flatnonzero(result[:,11]!=0)
    if len(failed) and not return_failures:
        raise RuntimeError(dict(operator='native exact Euler star',failed_records=failed[:16].tolist(),statuses=result[failed[:16],11].tolist(),fallback=False))
    return result
