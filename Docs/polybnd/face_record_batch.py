"""One outer evaluation's array preparation/native exact-star/ordered records."""
import time
import numpy as np
from native_euler_star import star_batch,provenance

def build(cache,edges,states,x,gen,frozen,centers,grad,boundary_velocity,gamma,alpha,expose_arrays=False):
    started=time.perf_counter();cpu=time.process_time();timers={}
    point=cache['points'];normal=cache['normals'];owners=cache['owners'];i=owners[:,0];j=owners[:,1];inside=j>=0
    count=len(point)
    if 'record_scale' in cache:
        # Saved-node oracle supplies the original measures directly. No new
        # quadrature or geometric inversion is performed in this audit path.
        scale=np.asarray(cache['record_scale']);distance=np.asarray(cache['record_distance'])
    else:
        sizes=np.diff(cache['offsets'])
        edge_lengths=np.array([e[5] for e in edges]);edge_distance=np.array([np.linalg.norm(x[e[1]]-x[e[0]]) if e[1]>=0 else 0. for e in edges])
        scale=np.concatenate([length*np.asarray(nodes)[:,1] for length,nodes in zip(edge_lengths,cache['edge_nodes'])])
        distance=np.repeat(edge_distance,sizes)
    wn=np.zeros(count)
    left_rho=cache['rho'][:,0];right_rho=cache['rho'][:,1].copy()
    left_mean=cache['means'][:,0];right_mean=cache['means'][:,1].copy()
    left_samples=cache['samples'][:,0];right_samples=cache['samples'][:,1].copy()
    left_xi=cache['coordinates'][:,0];right_xi=cache['coordinates'][:,1].copy()
    P=np.array([s['P'] for s in states])
    left_P=P[i]+(grad[i,None,:]@(point-centers[i])[:,:,None])[:,0,0]
    right_P=np.empty(count)
    ii=i[inside];jj=j[inside]
    right_P[inside]=P[jj]+(grad[jj,None,:]@(point[inside]-centers[jj])[:,:,None])[:,0,0]
    wn[inside]=(((x[jj]-x[ii])[:,None,:]@gen[jj,:,None])[:,0,0]-((gen[jj]-gen[ii])[:,None,:]@(point[inside]-x[ii])[:,:,None])[:,0,0])/distance[inside]
    exterior=~inside
    if exterior.any():
        external_velocity=np.array([boundary_velocity(p) for p in point[exterior]])
        wn[exterior]=(external_velocity[:,None,:]@normal[exterior,:,None])[:,0,0]
        dot=(left_mean[exterior,None,:]@normal[exterior,:,None])[:,0,0]
        right_mean[exterior]=left_mean[exterior]-2*(dot-wn[exterior])[:,None]*normal[exterior]
        sample_dot=(left_samples[exterior]@normal[exterior,:,None])[:,:,0]
        right_samples[exterior]=left_samples[exterior]-2*(sample_dot-wn[exterior,None])[:,:,None]*normal[exterior,None,:]
        right_rho[exterior]=left_rho[exterior];right_P[exterior]=left_P[exterior];right_xi[exterior]=alpha
    if not np.isfinite(left_P+right_P).all() or np.any(left_P<=0) or np.any(right_P<=0):raise RuntimeError('nonpositive/nonfinite reconstructed pressure')
    skip=frozen[i].copy();right_frozen=np.ones(count,dtype=bool);right_frozen[inside]=frozen[j[inside]];skip &= right_frozen
    ul=(left_mean[:,None,:]@normal[:,:,None])[:,0,0];ur=(right_mean[:,None,:]@normal[:,:,None])[:,0,0]
    inputs=np.column_stack([left_rho,ul,left_P,right_rho,ur,right_P,np.full(count,gamma)])
    timers['record_array_prepare']=dict(wall_seconds=time.perf_counter()-started,process_cpu_seconds=time.process_time()-cpu)
    started=time.perf_counter();cpu=time.process_time();out=np.full((count,13),np.nan);ps=np.zeros(count);sm=np.zeros(count)
    out[~skip]=star_batch(inputs[~skip]);ps[~skip]=out[~skip,0];sm[~skip]=out[~skip,1]
    timers['exact_star_native']=dict(wall_seconds=time.perf_counter()-started,process_cpu_seconds=time.process_time()-cpu)
    started=time.perf_counter();cpu=time.process_time()
    records=[(int(i[k]),int(j[k]),scale[k],wn[k],distance[k],left_samples[k],right_samples[k],float(left_rho[k]),float(right_rho[k]),normal[k],point[k],left_P[k],right_P[k],ps[k],sm[k],0.,left_xi[k],right_xi[k]) for k in range(count)]
    timers['record_tuple_pack']=dict(wall_seconds=time.perf_counter()-started,process_cpu_seconds=time.process_time()-cpu)
    audit=dict(inputs=inputs,star=out,skip=skip,pressure_error_bound=out[:,7],speed_error_bound=out[:,8])
    if expose_arrays:
        audit['transport_arrays']=dict(vl=left_samples,vr=right_samples,rl=left_rho,rr=right_rho,pl=left_P,pr=right_P,ps=ps,sm=sm,xil=left_xi,xir=right_xi,base=wn)
    return records,dict(nodes=count,star_calls=int((~skip).sum()),frozen_skipped=int(skip.sum()),
        timers=timers,provenance=provenance(),cache_scope='this RHS outer only'),audit
