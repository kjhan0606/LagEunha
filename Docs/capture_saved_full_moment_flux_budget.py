"""One saved t=.025 RHS with compact signed budgets; no time evolution.

Read-only instrumentation delegates the unchanged captured-node and GPU stage
operators. Original state, density fits, geometry and weight-rate seed are reused.
"""
import argparse,fcntl,hashlib,json,os,sys
from pathlib import Path
import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def arrays(path):
    with np.load(path) as f:return {k:f[k] for k in f.files}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    base=Path('/gpfs/kjhan/LagForce');fast=base/'Codex_full_disk_fast_modules_v1_20261019'
    cp=base/'Codex_full_disk_continuous_two_dt_to_t1_20261019/coarse/checkpoints/coarse_0005'
    old=cp.parent.parent/'steps/coarse_0006_start'
    files=[cp/'accepted.npz',cp/'geometry/mesh.npz',cp/'geometry/full_exact_supports.json',cp/'density/cells.jsonl',old/'report.json',old/'preparation.json']
    hashes={str(f):sha(f) for f in files};state=arrays(cp/'accepted.npz');mesh=arrays(cp/'geometry/mesh.npz')
    assert float(state['time'])==.025
    native=json.loads((cp/'geometry/full_exact_supports.json').read_text())
    fits=[json.loads(x) for x in (cp/'density/cells.jsonl').read_text().splitlines()]
    key=sha(cp/'accepted.npz')+':'+sha(cp/'geometry/full_exact_supports.json')
    prep=json.loads((old/'preparation.json').read_text());assert prep['input_state_sha256']==sha(cp/'accepted.npz')
    assert prep['geometry_support_sha256']==sha(cp/'geometry/full_exact_supports.json')
    # This diagnostic does not repeat geometry construction or density fitting.
    sys.path.insert(0,str(fast))
    import full_disk_stage_memory as stage
    import full_disk_velocity_memory as velocity
    from moment_compiled_backend import Backend
    original_capture=stage.capture;captured=[];n=len(state['ids'])
    def observed_capture(geometry,transport,values,states,gamma):
        result=original_capture(geometry,transport,values,states,gamma)
        if not captured:
            compact={k:np.zeros(n,np.longdouble) for k in ('Vdot','mass_rate','advected_internal_energy_rate','advected_kinetic_energy_rate','pressure_power_rate')}
            owners=result['owners'];node=result['node_inputs'];wn=values[2]
            for start in range(0,len(node),65536):
                end=min(start+65536,len(node));x=node[start:end].astype(np.longdouble);i,j=owners[start:end].T;inside=j>=0
                wl,wr=x[:,21:25],x[:,25:29];vl=x[:,5:13].reshape(-1,4,2);vr=x[:,13:21].reshape(-1,4,2)
                ds=x[:,4];fm=ds*(wl.sum(1)+wr.sum(1))
                streams={'Vdot':-ds*np.asarray(wn[start:end],np.longdouble),'mass_rate':fm,
                    'advected_internal_energy_rate':ds*(wl.sum(1)*x[:,29]+wr.sum(1)*x[:,30]),
                    'advected_kinetic_energy_rate':ds*((wl*.5*(vl*vl).sum(2)).sum(1)+(wr*.5*(vr*vr).sum(2)).sum(1)),
                    'pressure_power_rate':ds*x[:,31]*x[:,32]}
                for name,v in streams.items():
                    np.add.at(compact[name],i,-v);np.add.at(compact[name],j[inside],v[inside])
            np.savez(args.output/'signed_flux_volume.npz',ids=state['ids'],**compact)
            captured.append(dict(nodes=len(node),total_volume_rate=float(compact['Vdot'].sum())))
        return result
    stage.capture=observed_capture
    v=velocity.prepare(state,native,mesh,fits,key,workers=8)
    lock=base/'Codex_full_disk_continuous_two_dt_to_t1_20261019/GPU.lock'
    backend=Backend(base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so','cuda',1)
    class Locked:
        path=backend.path;kind=backend.kind
        def stage(self,data):
            with lock.open('a') as h:
                fcntl.flock(h.fileno(),fcntl.LOCK_EX)
                return backend.stage(data)
    try:
        result=stage.evaluate(state,native,mesh,fits,v.prepared,v.gravity_nodes,key,Locked(),args.output/'rhs',provider_budget=1)
        saved={k:value for k,value in result.saved.items() if k not in ('incoming_Q','incoming_offsets')}
        np.savez(args.output/'stage_compact.npz',**saved)
    finally:
        del v.gravity_nodes;v.owner.close()
    before=json.loads((old/'report.json').read_text());after=result.report
    same_nodes=all(before[k]==after[k] for k in ('node_inputs_sha256','node_extra_sha256'))
    assert all(sha(f)==h for f,h in hashes.items())
    report=dict(status='completed_single_saved_state_rhs_budget_capture',time=.025,
        physical_providers=1,accepted_time_steps=0,original_sources_unchanged=True,
        source_sha256=hashes,source_operator_sha256=sha(stage.__file__),source_velocity_sha256=sha(velocity.__file__),
        actual_backend_sha256=sha(backend.path),input_weight_rate_seed_sha256=hashlib.sha256(state['weight_rate_seed'].tobytes()).hexdigest(),
        original_node_hashes_identical=same_nodes,capture=captured[0],
        limits='One instantaneous RHS, not a trajectory or finite-step volume-map audit. Longdouble signed flux sums use different accumulation order and are diagnostic. No factor/time update, tolerance change or state repair.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()
