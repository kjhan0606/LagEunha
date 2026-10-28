"""Typed native state serialization; no numerical update or acceptance.

All cells remain active and constrained in this full-disk model. Long-double
potential is transported without conversion to binary64. The native driver
owns epoch/time, candidates and commitment; Python only reads/writes arrays.
"""
import hashlib,os,struct
from pathlib import Path
import numpy as np

MAGIC=0x31544154534d464c
VERSION=1
HEADER=struct.Struct('<8Qd7Q')
assert HEADER.size==128
FIELDS=('ids','mass','state','factors','Psi','weight_rate_seed')

def load(path):
    path=Path(path)
    if path.is_file():
        with np.load(path) as f:return {k:f[k] for k in f.files}
    data=(path/'header.bin').read_bytes();assert len(data)==HEADER.size
    magic,version,n,epoch,token,length,ldbytes,ldmant,time,*reserved=HEADER.unpack(data)
    assert magic==MAGIC and version==VERSION and n>0 and length==20*n+4
    assert ldbytes==np.dtype(np.longdouble).itemsize and ldmant==np.finfo(np.longdouble).nmant+1 and not any(reserved)
    assert np.isfinite(time) and time>=0
    shapes=dict(ids=(n,),mass=(n,),state=(20*n+4,),factors=(n,2,2),Psi=(n,),weight_rate_seed=(n,))
    result={}
    for key,shape in shapes.items():
        dtype=np.dtype('<i8') if key=='ids' else np.dtype(np.longdouble) if key=='Psi' else np.dtype('<f8')
        file=path/(key+'.bin');assert file.stat().st_size==int(np.prod(shape))*dtype.itemsize,(key,shape)
        result[key]=np.fromfile(file,dtype=dtype).reshape(shape)
    assert np.all(np.diff(result['ids'])>0) and (result['mass']>0).all()
    assert all(np.isfinite(value).all() for value in result.values())
    result.update(time=np.array(time),frozen=np.zeros(n,bool),psi_constrained_mask=np.ones(n,bool))
    return result

def dump(path,state,epoch,token):
    path=Path(path);path.mkdir(exist_ok=False);n=len(state['ids'])
    assert not state['frozen'].any() and state['psi_constrained_mask'].all()
    assert state['Psi'].dtype==np.dtype(np.longdouble)
    header=HEADER.pack(MAGIC,VERSION,n,epoch,token,20*n+4,np.dtype(np.longdouble).itemsize,
        np.finfo(np.longdouble).nmant+1,float(state['time']),*([0]*7))
    for key in FIELDS:
        value=np.ascontiguousarray(state[key]);assert value.dtype==(np.dtype('<i8') if key=='ids' else np.dtype(np.longdouble) if key=='Psi' else np.dtype('<f8'))
        with (path/(key+'.bin')).open('xb') as f:value.tofile(f)
    with (path/'header.bin').open('xb') as f:f.write(header)
    observed=load(path)
    assert all(np.array_equal(observed[k],state[k]) for k in observed), 'Serialization changed state'
    return path

def digest(path):
    path=Path(path);h=hashlib.sha256()
    for file in [path] if path.is_file() else [path/'header.bin',*[path/(key+'.bin') for key in FIELDS]]:
        with file.open('rb') as f:
            while data:=f.read(1024*1024):h.update(data)
    return h.hexdigest()
