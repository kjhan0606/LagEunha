"""Batch the existing shared face moment/covariance flux accumulation.

The scalar code constructs the same face states, HLLC fan, mass root and
weights. This module only groups their algebraic products across records.
"""
import numpy as np
from kinetic_factor_rhs import ALPHA


def backend_name():
    import os
    return 'batched' if os.environ.get('LAGFORCE_VECTOR_MOMENT_FLUX','0')=='1' else 'scalar'


def accumulate(records,values,states,mass,gamma):
    n=len(states);mass=np.asarray(mass);i=np.array([r[0] for r in records],int);j=np.array([r[1] for r in records],int);count=len(i)
    scale=np.array([r[2] for r in records]);vl=np.array([r[5] for r in records]);vr=np.array([r[6] for r in records]);rl=np.array([r[7] for r in records]);rr=np.array([r[8] for r in records]);normal=np.array([r[9] for r in records]);point=np.array([r[10] for r in records]);pl=np.array([r[11] for r in records]);pr=np.array([r[12] for r in records]);ps=np.array([r[13] for r in records]);sm=np.array([r[14] for r in records])
    wl=np.array([v[0] for v in values]);wr=np.array([v[1] for v in values]);wn=np.array([v[2] for v in values])
    xil=np.array([r[16] if len(r)>16 else ALPHA for r in records]);xir=np.array([r[17] if len(r)>17 else ALPHA for r in records])
    owners=np.column_stack([i,j]).ravel();valid=owners>=0
    def paired(left,right,shape):
        out=np.zeros((n,)+shape);terms=np.stack([left,right],axis=1).reshape((2*count,)+shape);np.add.at(out,owners[valid],terms[valid]);return out
    # Preserve the small per-record reduction cores of the scalar oracle.
    def sum4(w):return ((w[:,0]+w[:,1])+w[:,2])+w[:,3]
    fm=sum4(wl)+sum4(wr)
    fp=(wl[:,None,:]@vl)[:,0,:]+(wr[:,None,:]@vr)[:,0,:]+ps[:,None]*normal
    Pi=vl.transpose(0,2,1)@(wl[:,:,None]*vl)+vr.transpose(0,2,1)@(wr[:,:,None]*vr)
    el=.5*np.sum(vl*vl,axis=2)+pl[:,None]/((gamma-1)*rl[:,None]);er=.5*np.sum(vr*vr,axis=2)+pr[:,None]/((gamma-1)*rr[:,None])
    energy=(wl[:,None,:]@el[:,:,None])[:,0,0]+(wr[:,None,:]@er[:,:,None])[:,0,0]+ps*sm
    flux=scale[:,None]*np.column_stack([fm,point*fm[:,None],(point[:,:,None]*point[:,None,:]*fm[:,None,None]).reshape(-1,4),fp,(point[:,:,None]*fp[:,None,:]).reshape(-1,4),Pi.reshape(-1,4),energy])
    dr=paired(-flux,flux,(18,));area=paired(scale*wn,-scale*wn,())
    traction=scale[:,None]*ps[:,None]*normal;force=paired(-traction,traction,(2,))
    internal=j>=0
    center=np.array([s['c'] for s in states]);hleft=-(point-center[i])[:,:,None]*traction[:,None,:]
    hright=np.zeros_like(hleft)
    hright[internal]=(point[internal]-center[j[internal]])[:,:,None]*traction[internal,None,:]
    H=paired(hleft,hright,(2,2))
    kl=(xil.transpose(0,2,1)@(wl[:,:,None]*xil))*(scale/mass[i])[:,None,None]
    kr=np.zeros_like(kl)
    xr=xir[internal]
    kr[internal]=(xr.transpose(0,2,1)@((-wr[internal])[:,:,None]*xr))*(scale[internal]/mass[j[internal]])[:,None,None]
    K=paired(kl,kr,(2,2))
    velocity=np.array([s['v'] for s in states]);B=np.array([s['B'] for s in states])
    affine_left=velocity[i]+(B[i]@(point-center[i])[:,:,None])[:,:,0]
    qleft=np.sqrt(-scale[:,None]*wr)[:,None,:]*np.swapaxes(vr-affine_left[:,None,:],1,2)
    qright=np.zeros_like(qleft)
    affine_right=velocity[j[internal]]+(B[j[internal]]@(point[internal]-center[j[internal]])[:,:,None])[:,:,0]
    qright[internal]=np.sqrt(scale[internal,None]*wl[internal])[:,None,:]*np.swapaxes(vl[internal]-affine_right[:,None,:],1,2)
    blocks=np.stack([qleft,qright],axis=1).reshape(2*count,2,4)[valid];receivers=owners[valid];order=np.argsort(receivers,kind='stable');offset=np.r_[0,np.cumsum(np.bincount(receivers,minlength=n))]
    incoming=[]
    for a,b in zip(offset[:-1],offset[1:]):incoming.append([np.transpose(blocks[order[a:b]],(1,0,2)).reshape(2,-1)])
    boundary=j<0;ledger_terms=-scale[boundary,None]*np.column_stack([fp[boundary],energy[boundary],point[boundary,0]*fp[boundary,1]-point[boundary,1]*fp[boundary,0]])
    ledger=np.cumsum(ledger_terms,axis=0)[-1] if len(ledger_terms) else np.zeros(4)
    return dr,ledger,force,H,K,incoming,area
