"""Lossless, pickle-free storage of actual 3D midpoint paths and endpoints."""
from pathlib import Path
import numpy as np

KEYS=('t','dt','sites','weights','raw','factors','dx0','dw0','dr0','dxmid','dwmid','drmid')


def save(path,record,endpoint):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    d=record['diag_mid'];Q=d['cov_Qfactor'];offsets=np.r_[0,np.cumsum([q.shape[1] for q in Q])]
    arrays={k:np.asarray(record[k]) for k in KEYS}
    for key in ('face_integration','face_order','pressure_reconstruction','stage_face_orders'):
        if key in record:arrays[key]=np.asarray(record[key])
    arrays.update(cov_A=np.asarray(d['cov_A']),cov_K=np.asarray(d['cov_K']),cov_linearization_factors=np.asarray(d['cov_linearization_factors']),cov_Q=np.concatenate(Q,axis=1),offsets=offsets)
    for key,value in endpoint.items():arrays['end_'+key]=np.asarray(value)
    temp=Path(str(path)+'.tmp')
    with temp.open('wb') as f:np.savez_compressed(f,**arrays)
    temp.replace(path)


def load(path):
    with np.load(path,allow_pickle=False) as data:
        record={k:data[k].copy() for k in KEYS};offsets=data['offsets'];Q=data['cov_Q'];N=len(record['sites'])
        if len(offsets)!=N+1 or offsets[0]!=0 or offsets[-1]!=Q.shape[1] or np.any(np.diff(offsets)<0):raise ValueError('invalid covariance offsets')
        record['t']=float(record['t']);record['dt']=float(record['dt'])
        for key in ('face_integration','face_order','pressure_reconstruction'):
            if key in data.files:record[key]=data[key].item()
        if 'stage_face_orders' in data.files:record['stage_face_orders']=data['stage_face_orders'].tolist()
        record['diag_mid']={k:data[k].copy() for k in ('cov_A','cov_K','cov_linearization_factors')}
        record['diag_mid']['cov_Qfactor']=[Q[:,offsets[i]:offsets[i+1]].copy() for i in range(N)]
        endpoint={k[4:]:data[k].copy() for k in data.files if k.startswith('end_')}
    return record,endpoint
