"""Native exact non-arc area bisection; true-arc cells retain Python geometry."""
import ctypes
from pathlib import Path

class ExactPriorBisection:
    def __init__(self, library):
        self.library=Path(library).resolve(strict=True)
        self.lib=ctypes.CDLL(str(self.library))
        self.lib.exact_prior_bisect64.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_char_p),
            *([ctypes.c_double]*6),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_uint64)]
        self.lib.exact_prior_bisect64.restype=ctypes.c_int

    def bisect(self,cell,normal,lo,hi,target):
        if cell.arcs:raise ValueError('Native bisection requires an arc-free cell')
        strings=[str(v).encode() for point in cell.vertices for v in point]
        xy=(ctypes.c_char_p*len(strings))(*strings)
        result=ctypes.c_double();witness=ctypes.c_uint64()
        code=self.lib.exact_prior_bisect64(len(cell.vertices),xy,float(normal[0]),float(normal[1]),
            float(lo),float(hi),float(target),float(cell.radius),ctypes.byref(result),ctypes.byref(witness))
        if code:raise RuntimeError(dict(status='native_exact_prior_failure',code=code,owner=cell.owner))
        return result.value,witness.value
