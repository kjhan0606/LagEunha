"""Typed native kinetic-root primitives; sparse solve and audits stay unchanged."""
import ctypes as C
import hashlib
import os
from pathlib import Path

import numpy as np


class Prepared:
    def __init__(self, data):
        self.arrays = [np.ascontiguousarray(data[k], dtype=np.int64 if k in ('i', 'j') else np.float64)
                       for k in ('i', 'j', 'scale', 'base', 'distance', 'vl', 'vr', 'rl', 'rr')]
        self.n = len(self.arrays[0])
        assert all(x.shape == ((self.n, 4) if k in (5, 6) else (self.n,))
                   for k, x in enumerate(self.arrays))
        self.pointers = [x.ctypes.data_as(C.c_void_p) for x in self.arrays]
        i,j=self.arrays[:2];self.inside=j>=0
        a,b=i[self.inside],j[self.inside]
        self.jacobian_rows=np.r_[a,b,a,b]
        self.jacobian_columns=np.r_[a,b,b,a]


_library = None


def evaluate(prepared, count, z):
    global _library
    if _library is None:
        path = Path(os.environ['LAGFORCE_NATIVE_TRANSPORT_ROOT_LIB']).resolve()
        _library = C.CDLL(str(path))
        _library.lag_transport_root.argtypes = [C.c_int64, C.c_int64] + [C.c_void_p]*18
        _library.lag_transport_root.restype = C.c_int
    z = np.ascontiguousarray(z, dtype=np.float64)
    assert z.shape == (count,)
    n = prepared.n
    wl, wr = np.empty((n, 4)), np.empty((n, 4))
    wn, fm, primitive, rho = [np.empty(n) for _ in range(4)]
    residual = np.empty(count)
    weight = np.empty(n)
    outputs = [wl, wr, wn, fm, primitive, rho, residual, weight]
    status = _library.lag_transport_root(n, count, *prepared.pointers,
        z.ctypes.data_as(C.c_void_p), *(x.ctypes.data_as(C.c_void_p) for x in outputs))
    if status:
        raise RuntimeError(f'Native transport-root primitive failed: {status}')
    return dict(wl=wl, wr=wr, wn=wn, fm=fm, primitive=primitive,
                rho=rho, residual=residual, weight=weight)


def provenance():
    path = Path(os.environ['LAGFORCE_NATIVE_TRANSPORT_ROOT_LIB']).resolve()
    return dict(library=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                accumulation='Two canonical-order owner sums, then subtraction',
                sparse_solver_changed=False, tolerances_changed=False)
