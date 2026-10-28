"""Audit-only export of existing converged kinetic records and Python flux."""
import numpy as np

def capture(records,values,gamma,states=None):
    inputs=np.empty((len(records),33));expected=np.empty((len(records),18));owners=np.empty((len(records),2),dtype=np.int64)
    extra=np.zeros((len(records),20)) if states is not None else None
    for k,(rec,(wl,wr,wn)) in enumerate(zip(records,values)):
        i,j,scale,_,_,vl,vr,rl,rr,n,point,pl,pr,ps,sm,*_=rec
        eL=pl/((gamma-1)*rl);eR=pr/((gamma-1)*rr)
        inputs[k]=np.r_[point,n,scale,vl.ravel(),vr.ravel(),wl,wr,eL,eR,ps,sm]
        fm=sum(wl)+sum(wr);fp=wl@vl+wr@vr+ps*n
        Pi=vl.T@(wl[:,None]*vl)+vr.T@(wr[:,None]*vr)
        fe=wl@(.5*np.sum(vl*vl,axis=1)+eL)+wr@(.5*np.sum(vr*vr,axis=1)+eR)+ps*sm
        expected[k]=scale*np.r_[fm,point*fm,(np.outer(point,point)*fm).ravel(),fp,np.outer(point,fp).ravel(),Pi.ravel(),fe]
        owners[k]=(i,j)
        if extra is not None:
            if len(rec)!=18:raise ValueError('quadratic coordinates required for stage capture')
            extra[k,:8]=np.asarray(rec[16]).ravel();extra[k,8:16]=np.asarray(rec[17]).ravel()
            for side,receiver in enumerate((i,j)):
                if receiver>=0:
                    s=states[receiver]
                    extra[k,16+2*side:18+2*side]=s['v']+s['B']@(point-s['c'])
    result=dict(inputs=inputs,expected=expected,owners=owners)
    if extra is not None:result['extra']=extra
    return result
