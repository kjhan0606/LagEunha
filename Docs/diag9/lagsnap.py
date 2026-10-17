"""khout.*.dat reader. ps=576: VoroQ float (before 2026-10-02); ps=624: VoroQ double."""
import struct, numpy as np
_Q = ['den','pressure','avgP','ie','die','csound','dt','volume','w2','w2hydro','w2old','w2ceil']
def read(path):
    with open(path,'rb') as f:
        _,n,ps,_ = struct.unpack('iiii', f.read(16))
        a = np.frombuffer(f.read(n*ps), dtype=np.uint8).reshape(n, ps)
        t = struct.unpack('d', f.read()[:8])[0]
    if ps == 576: qt, qs = np.float32, 4
    elif ps == 624: qt, qs = np.float64, 8
    else: raise ValueError(f'unknown particle size {ps}')
    g = lambda o, dt=np.float64, s=8: a[:, o:o+s].copy().view(dt).ravel()
    out = dict(t=t, n=n, indx=g(8, np.int64)); out['id'] = out['indx'] & ((1 << 56) - 1)   # PINDX
    for k, nm in enumerate(['x','y','z','mass','vx','vy','vz','ax','ay','az']): out[nm] = g(24 + 8*k)
    for k, nm in enumerate(_Q): out[nm] = g(104 + qs*k, qt, qs).astype(np.float64)
    return out
