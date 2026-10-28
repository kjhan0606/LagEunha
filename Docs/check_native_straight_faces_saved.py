"""Initial/moving true-disk face component oracles; no massroot or fluid RHS."""
import argparse
import hashlib
import importlib.util
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
    args=parser.parse_args();args.output.mkdir(exist_ok=False)
    own=Path(__file__).resolve().parent;base=Path('/gpfs/kjhan/LagForce')
    sys.path.insert(0,str(base/'Codex_full_disk_fast_modules_v1_20261019'))
    import full_disk_velocity_memory as velocity
    velocity.model.GAMMA=float(np.float32(1.666667))
    import curved_face_nodes_fast as original
    from same_stage_face_cache import SameStageFaceCache
    from curved_disk_cell import DiskCell
    from compact_potential_density import CompactPotentialDensity,VolumeCache
    spec=importlib.util.spec_from_file_location('native_geometry_component',own/'curved_face_nodes_native.py')
    candidate=importlib.util.module_from_spec(spec);spec.loader.exec_module(candidate)
    import native_straight_faces
    assert os.environ['OMP_NUM_THREADS']==os.environ['OPENBLAS_NUM_THREADS']=='1'
    reports={};passed=True
    for label in ('initial','moving'):
        if label=='initial':
            state_path=base/'Codex_full_curved_kepler_IC_20261019/initial.npz'
            mesh_path=base/'Codex_full_curved_disk_mesh_20261019/mesh.npz'
            support_path=base/'Codex_circular_sparse_support_full_20261019/full_exact_supports.json'
            fit_path=base/'Codex_full_disk_initial_density_20261019/cells.jsonl'
        else:
            cp=base/'Codex_full_disk_continuous_two_dt_to_t1_20261019/coarse/checkpoints/coarse_0005'
            state_path=cp/'accepted.npz';mesh_path=cp/'geometry/mesh.npz'
            support_path=cp/'geometry/full_exact_supports.json';fit_path=cp/'density/cells.jsonl'
        source={str(p):sha(p) for p in (state_path,mesh_path,support_path,fit_path)}
        state=read(state_path);mesh=read(mesh_path);support=json.loads(support_path.read_text())
        fits=[json.loads(line) for line in fit_path.read_text().splitlines()]
        assert len(fits)==len(state['ids']) and all(row['status']=='passed' for row in fits)
        stage_key=sha(state_path)+':'+sha(support_path)
        print(json.dumps(dict(label=label,status='preparing_saved_velocity_metadata',physical_providers=0)),flush=True)
        if label=='initial':
            p=base/'Codex_full_disk_initial_velocity_gravity_20261019/prepared.npz'
            source[str(p)]=sha(p);prepared=read(p)
        else:
            v=velocity.prepare(state,support,mesh,fits,stage_key,workers=8)
            prepared=v.prepared;del v.gravity_nodes;v.owner.close()
        assert all(np.array_equal(prepared[key],state[key]) for key in ('state','mass','factors','Psi','ids'))
        mass=state['mass'];ids=state['ids'];L=state['factors'];a=state['state'][:-4].reshape(-1,20)
        raw=np.column_stack([mass,a[:,3:]]);positions=a[:,:2]
        cache=VolumeCache(1);cache.begin_stage(stage_key);cells=[];states=[]
        for i,row in enumerate(fits):
            cell=DiskCell.from_native(support['cells'][i],2.5)
            density=CompactPotentialDensity(cell,raw[i],row['prior'],row['parameter'],velocity.potential,cache,stage_key,int(ids[i]))
            fluid=velocity.model.decode(cell,raw[i],density=density);fluid['R']=L[i]@L[i].T
            fluid['velocity_quadratic']={key:prepared[key][i] for key in ('center','scale','projection','coefficient','standard_coefficient','residual_map')}
            cells.append(cell);states.append(fluid)
        stage_cache=SameStageFaceCache(stage_key)
        arguments=(cells,ids,positions,mass,L,states,prepared['generator'],state['weight_rate_seed'],[r['prior'] for r in fits])
        begin=time.perf_counter();expected=original.build_geometry(*arguments,stage_cache=stage_cache);python_seconds=time.perf_counter()-begin
        begin=time.perf_counter();actual=candidate.build_geometry(*arguments,stage_cache=stage_cache);native_seconds=time.perf_counter()-begin
        comparisons={}
        for key in ('points','normals','owners','record_scale','record_distance','offsets'):
            x,y=expected[key],actual[key];okay=bool(np.array_equal(x,y));passed &= okay
            comparisons[key]=dict(exact_value=okay,reference_shape=list(x.shape),actual_shape=list(y.shape),
                maximum_difference=float(abs(x-y).max(initial=0.)) if x.shape==y.shape else None)
            if not okay and x.shape==y.shape:
                idx=np.flatnonzero(x.ravel()!=y.ravel())[:64]
                np.savez(args.output/(label+'_'+key+'_failure.npz'),indices=idx,reference=x.ravel()[idx],actual=y.ravel()[idx])
        for key in ('keys','arc_root_audit','geometry_stats'):
            okay=expected[key]==actual[key];comparisons[key]=dict(exact_value=okay);passed &= okay
        assert all(sha(Path(p))==digest for p,digest in source.items())
        reports[label]=dict(cells=len(ids),source_sha256=source,comparisons=comparisons,
            python_first_build_seconds=python_seconds,native_build_using_existing_stage_cache_seconds=native_seconds,
            timing_limit='First build versus same-stage cached build is NOT an apples-to-apples speedup measurement')
        print(json.dumps(dict(label=label,passed=passed,report=reports[label])),flush=True)
        if not passed:break
        del expected,actual,prepared,cells,states,stage_cache
    report=dict(status='passed_saved_straight_face_component_oracles' if passed else 'failed_saved_straight_face_component_oracles',
        stages=reports,provenance=native_straight_faces.provenance(),physical_providers=0,accepted_time_steps=0,
        density_Newton_fits=0,moving_velocity_preparations=1 if 'moving' in reports else 0,
        limits='True-disk saved-state component reconstruction only; full provider replay and fair exclusive timing remain required.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    assert passed,'Preserved first failure; no tolerance change or fallback'


if __name__=='__main__':main()
