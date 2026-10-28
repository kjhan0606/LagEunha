"""Two saved initial/moving transport component oracles, zero physical RHS."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import native_transport_root as native


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(exist_ok=False)
    root=Path('/gpfs/kjhan/LagForce');reports={};passed=True
    for label in ('stage0','midpoint'):
        source=root/'Codex_compiled_midpoint_once_20261018'/label/'inputs.npz'
        density_source=root/'Codex_face_record_native_saved_20261019'/(label+'.npz')
        before={str(p):sha(p) for p in (source,density_source)}
        with np.load(source) as f:
            nodes=f['node_inputs'];owners=f['owners'];state=f['state'];gen=f['generator']
            z=f['weight_rates'];count=len(f['mass'])
        with np.load(density_source) as f:
            riemann=f['inputs'];assert np.array_equal(f['owners'],owners)
        i,j=owners.T;inside=j>=0;point=nodes[:,:2];normal=nodes[:,2:4]
        position=state[:-4].reshape(count,20)[:,:2];distance=np.zeros(len(i));base=np.zeros(len(i))
        ii,jj=i[inside],j[inside];dx=position[jj]-position[ii]
        distance[inside]=np.linalg.norm(dx,axis=1)
        base[inside]=((dx[:,None,:]@gen[jj,:,None])[:,0,0]-
            ((gen[jj]-gen[ii])[:,None,:]@(point[inside]-position[ii])[:,:,None])[:,0,0])/distance[inside]
        data=dict(i=i,j=j,scale=nodes[:,4],base=base,distance=distance,
            vl=(nodes[:,5:13].reshape(-1,4,2)@normal[:,:,None])[:,:,0],
            vr=(nodes[:,13:21].reshape(-1,4,2)@normal[:,:,None])[:,:,0],
            rl=riemann[:,0],rr=riemann[:,3])
        begin=time.perf_counter()
        wn=base.copy();wn[inside]+=(z[i[inside]]-z[j[inside]])/(2*distance[inside])
        sl=data['vl']-wn[:,None];sr=data['vr']-wn[:,None]
        wl=.25*data['rl'][:,None]*np.maximum(sl,0.)
        wr=.25*data['rr'][:,None]*np.minimum(sr,0.)
        fm=data['scale']*(wl.sum(1)+wr.sum(1))
        primitive=data['scale']*(abs(wl).sum(1)+abs(wr).sum(1))
        residual=np.bincount(i,weights=fm,minlength=count)-np.bincount(j[inside],weights=fm[inside],minlength=count)
        rho=.25*(data['rl']*np.where(sl>0,1.,np.where(sl<0,0.,.5)).sum(1)+data['rr']*np.where(sr<0,1.,np.where(sr>0,0.,.5)).sum(1))
        weight=np.zeros(len(i));weight[inside]=data['scale'][inside]*rho[inside]/(2*distance[inside])
        reference=dict(wl=wl,wr=wr,wn=wn,fm=fm,primitive=primitive,residual=residual,rho=rho,weight=weight)
        python_seconds=time.perf_counter()-begin
        begin=time.perf_counter();prepared=native.Prepared(data);prepare_seconds=time.perf_counter()-begin
        begin=time.perf_counter();actual=native.evaluate(prepared,count,z);native_seconds=time.perf_counter()-begin
        comparisons={}
        for key,expected in reference.items():
            value=actual[key];okay=bool(np.array_equal(value,expected));passed &= okay
            comparisons[key]=dict(exact_value=okay,maximum_difference=float(abs(value-expected).max(initial=0.)))
            if not okay:
                idx=np.flatnonzero(value.ravel()!=expected.ravel())[:64]
                np.savez(args.output/(label+'_'+key+'_failure.npz'),indices=idx,expected=expected.ravel()[idx],actual=value.ravel()[idx])
        assert all(sha(Path(p))==digest for p,digest in before.items())
        reports[label]=dict(nodes=len(i),cells=count,input_sha256=before,comparisons=comparisons,
            python_primitive_seconds=python_seconds,native_prepare_seconds=prepare_seconds,native_primitive_seconds=native_seconds,
            scope='Component arithmetic on saved transport input; not a closed-root physical solve or provider comparison')
        if not passed:break
    report=dict(status='passed_saved_transport_component_oracles' if passed else 'failed_saved_transport_component_oracles',
        stages=reports,provenance=native.provenance(),physical_providers=0,accepted_time_steps=0,
        limits='No tolerance enlarged. Component timing excludes sparse solve and physics preparation; whole-provider replay still required.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    assert passed, 'First failed component evidence preserved; no fallback'


if __name__=='__main__':main()
