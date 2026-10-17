"""Continue the accepted local correction, preserving its original prefix.

Each accepted chunk is saved before starting the next. No old downstream
coefficients are reused after the corrected local endpoint.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from local_weight_defect import Problem, serialize, audit_window
from residual_controlled_path import integrate
from replay_dense_weight_path import segment_value
from dense_path_derivative import segment_derivative


def unpack(data):
    return [SimpleNamespace(t_old=s['t_old'],h=s['h'],
        y_old=np.asarray(s['y_old']),F=np.asarray(s['F'])) for s in data['segments']]


def continuity(segments, n):
    state_jump=0.;derivative_jump=0.
    for left,right in zip(segments,segments[1:]):
        np.testing.assert_allclose(left.t_old+left.h,right.t_old,rtol=0,atol=2e-15)
        state_jump=max(state_jump,float(np.max(abs(segment_value(left,right.t_old)-right.y_old))))
        derivative_jump=max(derivative_jump,float(np.max(abs(
            segment_derivative(left,right.t_old)[:n]-segment_derivative(right,right.t_old)[:n]))))
    return dict(max_state_jump=state_jump,max_weight_derivative_jump=derivative_jump)


def run(source, local_source):
    started=time.monotonic()
    source=Path(source);local_source=Path(local_source)
    saved=json.loads(source.with_name(source.name.replace('.dense.json','.path.json')).read_text())
    original=json.loads(source.read_text());local=json.loads(local_source.read_text())
    local_audit=json.loads(local_source.with_name(local_source.name.replace('.dense.json','.audit.json')).read_text())['result']
    if not local_audit['gate'] or local_audit['target']!=1e-9:
        raise RuntimeError('local predecessor gate failed')
    old=unpack(original);corrected=unpack(local)
    a,b=local['window']
    prefix=[s for s in old if s.t_old+s.h<=a+2e-15]
    if not prefix or abs(prefix[-1].t_old+prefix[-1].h-a)>2e-15:
        raise RuntimeError('prefix does not meet local start')
    np.testing.assert_allclose(segment_value(prefix[-1],a),local['initial'],rtol=0,atol=2e-15)
    np.testing.assert_array_equal(corrected[0].y_old,local['initial'])
    np.testing.assert_array_equal(segment_value(corrected[-1],b),local['end'])
    problem=Problem(saved)
    segments=prefix+corrected
    y=np.asarray(local['end']);s=b
    provenance=dict(original=str(source),local=str(local_source),
        original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        local_sha256=hashlib.sha256(local_source.read_bytes()).hexdigest(),
        prefix_segments=len(prefix),local_segments=len(corrected),resume_s=b,
        controller_target=1e-9,global_gate=1e-8,chunk_max_width=1/32,
        note='original prefix retained; corrected local endpoint feeds new suffix')
    Path('continuation_provenance.json').write_text(json.dumps(provenance,indent=2))
    checkpoints=Path('checkpoints');checkpoints.mkdir(exist_ok=False)
    index=0;rejections=0
    # Initial checkpoint is enough to resume even before the first new chunk.
    def save_checkpoint():
        out=dict(s=s,end=y.tolist(),segments=[serialize(q) for q in segments],
                 ts=[q.t_old for q in segments]+[s],rejections=rejections,
                 chunk=index,provenance=provenance)
        with (checkpoints/f'chunk_{index:04d}.json').open('x') as f:json.dump(out,f)
    save_checkpoint()
    while s<1.:
        right=min(1.,s+1/32)
        def progress(row):
            with open('attempts.jsonl','a') as f:f.write(json.dumps(dict(chunk=index+1,**row))+'\n')
        new,end,attempts=integrate(problem.fun,y,s,right,problem.check,1e-9,progress=progress)
        np.testing.assert_array_equal(new[0].y_old,y)
        segments.extend(new);y=end;s=right;index+=1
        rejections+=sum(not row['accepted'] for row in attempts)
        save_checkpoint()
        print(json.dumps(dict(stage='checkpoint',chunk=index,s=s,segments=len(segments),
            rejections=rejections,seconds=time.monotonic()-started)),flush=True)
    joins=continuity(segments,problem.p.N)
    Path('continued.dense.json').write_text(json.dumps(dict(ts=[q.t_old for q in segments]+[1.],
        segments=[serialize(q) for q in segments],end=y.tolist(),provenance=provenance)))
    print(json.dumps(dict(stage='integration_complete',**joins)),flush=True)
    result,rows=audit_window(problem,segments)
    prefix_count=len(prefix)
    result.update(joins)
    result.update(status='completed',seconds=time.monotonic()-started,calls=problem.calls,
        chunks=index,rejections=rejections,resume_s=b,final_s=s,
        prefix_max_residual=max(r['residual'] for r in rows if r['segment']<prefix_count),
        corrected_max_residual=max(r['residual'] for r in rows if r['segment']>=prefix_count),
        ode_accumulated_relative_mass=float(np.max(abs(y[problem.p.N:-4].reshape(problem.p.N,4)[:,0])/problem.p.mass)),
        limitation='frozen predictor, sampled path validation; not coupled hydro evolution')
    result['gate']=bool(result['gate'] and joins['max_state_jump']<1e-12
                        and joins['max_weight_derivative_jump']<1e-10)
    Path('continued.audit.json').write_text(json.dumps(dict(result=result,samples=rows),indent=2))
    Path('result.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source');parser.add_argument('local_source')
    args=parser.parse_args();run(args.source,args.local_source)
