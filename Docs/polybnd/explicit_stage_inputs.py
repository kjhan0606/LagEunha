"""Existing compiled-stage ABI from actual explicit nodes, without replay flux."""
import numpy as np

def capture(geometry,transport,values,states,gamma):
    wl,wr,wn=values;p=geometry['points'];n=geometry['normals'];owners=geometry['owners'];count=len(p)
    v=transport;inputs=np.column_stack([p,n,geometry['record_scale'],v['vl'].reshape(count,8),v['vr'].reshape(count,8),wl,wr,v['pl']/((gamma-1)*v['rl']),v['pr']/((gamma-1)*v['rr']),v['ps'],v['sm']])
    extra=np.zeros((count,20));extra[:,:8]=v['xil'].reshape(count,8);extra[:,8:16]=v['xir'].reshape(count,8)
    center=np.array([s['c'] for s in states]);velocity=np.array([s['v'] for s in states]);B=np.array([s['B'] for s in states]);i=owners[:,0];j=owners[:,1];inside=j>=0
    extra[:,16:18]=velocity[i]+(B[i]@(p-center[i])[:,:,None])[:,:,0]
    extra[inside,18:20]=velocity[j[inside]]+(B[j[inside]]@(p[inside]-center[j[inside]])[:,:,None])[:,:,0]
    flat=owners.ravel();valid=flat>=0;receivers=flat[valid];order=np.argsort(receivers,kind='stable');codes=np.column_stack([-np.arange(1,count+1),np.arange(1,count+1)]).ravel()[valid]
    incidence=np.ascontiguousarray(codes[order],dtype=np.int64);offsets=np.r_[0,np.cumsum(np.bincount(receivers,minlength=len(states)))].astype(np.int64)
    return dict(node_inputs=np.ascontiguousarray(inputs),extra=extra,owners=owners,offsets=offsets,incidence=incidence)


def primitive_energy_momentum(inputs):
    vl=inputs[:,5:13].reshape(-1,4,2);vr=inputs[:,13:21].reshape(-1,4,2);wl=inputs[:,21:25];wr=inputs[:,25:29];ds=inputs[:,4]
    def sum4(w):return ((w[:,0]+w[:,1])+w[:,2])+w[:,3]
    fm=ds*(sum4(wl)+sum4(wr));el=.5*np.sum(vl*vl,axis=2)+inputs[:,29,None];er=.5*np.sum(vr*vr,axis=2)+inputs[:,30,None]
    energy=ds*((wl[:,None,:]@el[:,:,None])[:,0,0]+(wr[:,None,:]@er[:,:,None])[:,0,0]+inputs[:,31]*inputs[:,32])
    fp=(wl[:,None,:]@vl)[:,0,:]+(wr[:,None,:]@vr)[:,0,:]+inputs[:,31,None]*inputs[:,2:4]
    absolute_fp=ds[:,None]*((abs(wl)[:,None,:]@abs(vl))[:,0,:]+(abs(wr)[:,None,:]@abs(vr))[:,0,:]+abs(inputs[:,31,None])*abs(inputs[:,2:4]))
    absolute_energy=ds*(np.sum(abs(wl*el),axis=1)+np.sum(abs(wr*er),axis=1)+abs(inputs[:,31]*inputs[:,32]))
    absolute_angular=abs(inputs[:,0])*absolute_fp[:,1]+abs(inputs[:,1])*absolute_fp[:,0]
    return fm,energy,ds[:,None]*fp,np.r_[absolute_fp.sum(0),absolute_energy.sum(),absolute_angular.sum()]
