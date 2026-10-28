"""One baseline and one native saved moving provider; no time evolution."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):
    with np.load(path) as f:return {key:f[key] for key in f.files}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    b=Path('/gpfs/kjhan/LagForce');own=Path(__file__).resolve().parent
    component_paths=[b/'Codex_native_transport_root_component_oracle_v1_20261019/report.json',
        b/'Codex_native_straight_faces_component_oracle_v1_20261019/report.json']
    statuses=['passed_saved_transport_component_oracles','passed_saved_straight_face_component_oracles']
    for path,status in zip(component_paths,statuses):
        row=json.loads(path.read_text());assert row['status']==status
        assert row['physical_providers']==row['accepted_time_steps']==0
        assert len(row['stages'])==2
    face=json.loads(component_paths[1].read_text())
    assert all(row['cells']==19945 for row in face['stages'].values())
    for path,variable in zip(component_paths,('LAGFORCE_NATIVE_TRANSPORT_ROOT_LIB','LAGFORCE_NATIVE_STRAIGHT_FACE_LIB')):
        row=json.loads(path.read_text());assert sha(Path(os.environ[variable]))==row['provenance']['sha256']
    args.output.mkdir(exist_ok=False)
    sys.path.insert(0,str(b/'Codex_full_disk_fast_modules_v1_20261019'))
    import full_disk_stage_memory as baseline
    import full_disk_velocity_memory as velocity
    sys.path.insert(0,str(own))
    import native_full_disk_stage as candidate
    from moment_compiled_backend import Backend
    assert os.environ['OMP_NUM_THREADS']==os.environ['OPENBLAS_NUM_THREADS']=='1'
    cp=b/'Codex_full_disk_continuous_two_dt_to_t1_20261019/coarse/checkpoints/coarse_0005'
    sources=[cp/'accepted.npz',cp/'geometry/mesh.npz',cp/'geometry/full_exact_supports.json',cp/'density/cells.jsonl']
    hashes={str(p):sha(p) for p in sources};state=read(sources[0]);mesh=read(sources[1])
    support=json.loads(sources[2].read_text());fits=[json.loads(line) for line in sources[3].read_text().splitlines()]
    assert float(state['time'])==.025;key=hashes[str(sources[0])]+':'+hashes[str(sources[2])]
    v=velocity.prepare(state,support,mesh,fits,key,workers=8)
    backend=Backend(b/'Codex_compiled_midpoint_build_20261018/liblagmoment.so','cuda',1)
    lock=b/'Codex_full_disk_continuous_two_dt_to_t1_20261019/GPU.lock'
    class Locked:
        path=backend.path;kind=backend.kind
        def stage(self,data):
            with lock.open('a') as f:
                fcntl.flock(f.fileno(),fcntl.LOCK_EX)
                return backend.stage(data)
    manifest=dict(status='started_two_saved_provider_performance_replays',input_sha256=hashes,
        component_evidence_sha256={str(p):sha(p) for p in component_paths},
        accepted_time_steps=0,maximum_physical_providers=2,source_sha256={
            str(p):sha(p) for p in own.iterdir() if p.is_file() and p.suffix in ('.py','.cpp','.so')},
        actual_backend_sha256=sha(backend.path),OMP_threads=1,BLAS_threads=1,velocity_workers=8)
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    results={};elapsed={}
    try:
        for label,stage in [('baseline',baseline),('candidate',candidate)]:
            print(json.dumps(dict(label=label,status='physical_provider_start',time=.025)),flush=True)
            begin=time.perf_counter()
            result=stage.evaluate(state,support,mesh,fits,v.prepared,v.gravity_nodes,key,Locked(),args.output/label,provider_budget=2)
            elapsed[label]=time.perf_counter()-begin;results[label]=result
            print(json.dumps(dict(label=label,status='physical_provider_completed',seconds=elapsed[label])),flush=True)
        a,c=results['baseline'],results['candidate'];comparisons={};passed=True
        assert a.saved.keys()==c.saved.keys()
        for name,x in a.saved.items():
            y=c.saved[name];okay=bool(np.array_equal(x,y));passed &= okay
            comparisons[name]=dict(exact_value=okay,shape=list(np.shape(x)))
            if not okay:
                indices=np.flatnonzero(np.asarray(x).ravel()!=np.asarray(y).ravel())[:64]
                np.savez(args.output/(name+'_failure.npz'),indices=indices,baseline=np.asarray(x).ravel()[indices],candidate=np.asarray(y).ravel()[indices])
        for name in ('node_inputs_sha256','node_extra_sha256'):
            okay=a.report[name]==c.report[name];comparisons[name]=dict(exact_value=okay);passed &= okay
        history_a=[(r['nodes'],r['closed_root'],r['mass_residual'],r['weight_rate_change']) for r in a.report['outer']]
        history_c=[(r['nodes'],r['closed_root'],r['mass_residual'],r['weight_rate_change']) for r in c.report['outer']]
        same_history=history_a==history_c;passed &= same_history
        profile={label:dict(total_seconds=elapsed[label],geometry_seconds=sum(r['geometry_seconds'] for r in result.report['outer']),
            massroot_seconds=sum(r['massroot_seconds'] for r in result.report['outer']),peak_RSS_KiB=result.report['peak_RSS_KiB']) for label,result in results.items()}
        combined=lambda p:p['geometry_seconds']+p['massroot_seconds']
        speed=combined(profile['baseline'])/combined(profile['candidate'])
        assert all(sha(Path(p))==digest for p,digest in hashes.items())
        report=dict(status='passed_native_face_provider_equivalence' if passed else 'failed_native_face_provider_equivalence',
            physical_providers=2,accepted_time_steps=0,time=.025,cells=len(state['ids']),comparisons=comparisons,
            closed_root_history_exact_value=same_history,profile=profile,geometry_root_speed_ratio=speed,
            geometry_root_2x_goal_met=speed>=2,whole_provider_speed_ratio=elapsed['baseline']/elapsed['candidate'],
            source_inputs_unchanged=True,native_provenance=c.native_provenance,
            limits='One same-host saved-state pair. No time update, integrator convergence or long-time stability claim; speed goal failure is retained. RSS is cumulative process high-water, including retained baseline arrays, not independent per-provider peak.')
        (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
        assert passed,'First full-provider failure retained; no fallback or enlarged bounds'
    finally:
        del v.gravity_nodes;v.owner.close()


if __name__=='__main__':main()
