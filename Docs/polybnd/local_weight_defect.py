"""Diagnose and refine only the saved worst interval and its neighbors.

No splice into the old global trajectory: changing a window's endpoint would
require continuing the downstream solution. This local experiment tests the
controller before any such continuation is authorized by its numerical gate.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.optimize import root
from predicted_spacetime_flux import FluxAudit
from kinetic_factor_rhs import ALPHA
import moment_ale_evolve as base
from dense_path_derivative import segment_derivative
from replay_dense_weight_path import segment_value
from residual_controlled_path import integrate


def serialize(segment):
    return dict(t_old=float(segment.t_old), h=float(segment.h),
                y_old=segment.y_old.tolist(), F=segment.F.tolist())


class Problem:
    def __init__(self, saved):
        r = saved['result']
        self.audit = FluxAudit(r['kind'], r['n'], r['duration'])
        self.p = self.audit.path
        for actual, expected in zip(self.p.x, saved['generator_polynomials']):
            for poly, coefficients in zip(actual, expected):
                np.testing.assert_allclose(poly.coef, coefficients, rtol=0, atol=1e-14)
        self.p.x = [[P(c) for c in row] for row in saved['generator_polynomials']]
        self.calls = 0

    def weights(self, w):
        self.p.w = [P([float(v)]) for v in w]

    def constrained(self, t, w):
        self.weights(w)
        cache = {}
        def evaluate(z):
            key = z.tobytes()
            if cache.get('key') != key:
                self.calls += 1
                f, j = self.audit.rate(t, np.r_[z, 0.], return_jacobian=True)
                cache.update(key=key, f=f[:-4].reshape(self.p.N, 4)[:-1, 0], j=j[:-1, :-1])
            return cache
        answer = root(lambda z: evaluate(z)['f'], np.zeros(self.p.N-1),
                      jac=lambda z: evaluate(z)['j'], method='hybr', options={'xtol':1e-9})
        z = np.r_[answer.x, 0.]
        f = self.audit.rate(t, z)
        if not np.all(np.isfinite(z)) or not np.all(np.isfinite(f)) or np.max(abs(f[:-4].reshape(self.p.N,4)[:,0])) > 2e-14:
            raise RuntimeError('constraint solve residual failed')
        return z, f

    def fun(self, t, y):
        z, f = self.constrained(t, y[:self.p.N])
        return np.r_[z, f]

    def flux(self, segment, t):
        self.weights(segment_value(segment, t)[:self.p.N])
        flux = self.audit.rate(t, segment_derivative(segment, t)[:self.p.N])
        if not np.all(np.isfinite(flux)):
            raise RuntimeError('nonfinite path flux')
        return flux

    def residual(self, segment, t):
        return float(np.max(abs(self.flux(segment,t)[:-4].reshape(self.p.N,4)[:,0])/self.p.mass))

    def check(self, segment):
        return max(self.residual(segment,t) for t in np.linspace(segment.t_old, segment.t_old+segment.h,17))

    def branches(self, t, w, z):
        """Internal face speed signs: zero crossings switch kinetic upwinding."""
        self.weights(w)
        p = self.p
        x, weights = p.state(t)
        _, edges = base.geometry(x, weights, t*p.h)
        gx = (np.polynomial.legendre.leggauss(4)[0]+1)/2
        faces = []
        speeds = {}
        for i,j,a,b,n,length in edges:
            if j < 0:
                continue
            faces.append([int(i),int(j),float(length)])
            row = p.row(i,j)
            normal = np.array([row[0](t),row[1](t)])
            for q, g in enumerate(gx):
                point = a+g*(b-a)
                wn = (row[2].deriv()(t)-point@np.array([row[0].deriv()(t),row[1].deriv()(t)]))/np.linalg.norm(normal)/p.h
                wn += (z[i]-z[j])/(2*np.linalg.norm(x[j]-x[i])*p.h)
                for side, cell in enumerate((i,j)):
                    velocity = base.velocity(self.audit.states[cell],point)
                    samples = velocity+ALPHA@p.mid_factors[cell].T/np.sqrt(p.mass[cell])
                    for sample, speed in enumerate(samples@n-wn):
                        speeds[f'{i}:{j}:{q}:{side}:{sample}'] = float(speed)
        return faces, speeds


def audit_window(problem, segments):
    # 129 points include 65 validation nodes and their staggered midpoints.
    rows = []
    integrals = {q:np.zeros(4*problem.p.N+4) for q in (16,32)}
    for k,seg in enumerate(segments):
        for t in np.linspace(seg.t_old,seg.t_old+seg.h,129):
            rows.append(dict(segment=k,s=float(t),residual=problem.residual(seg,t)))
        for order in integrals:
            nodes,weights = np.polynomial.legendre.leggauss(order)
            for x,w in zip(nodes,weights):
                integrals[order] += seg.h*w/2*problem.flux(seg,seg.t_old+seg.h*(x+1)/2)
    p=problem.p
    mass=lambda f:float(np.max(abs(f[:-4].reshape(p.N,4)[:,0])/p.mass))
    worst=max(rows,key=lambda row:row['residual'])
    f=integrals[32]
    result=dict(worst=worst,samples=len(rows),segments=len(segments),
                mass_defect=mass(f),quadrature_mass_difference=mass(f-integrals[16]),
                ledger_error=float(np.max(abs(f[:-4].reshape(p.N,4).sum(0)-f[-4:]))))
    result['gate']=bool(worst['residual']<1e-8 and result['mass_defect']<1e-10
                        and result['quadrature_mass_difference']<1e-11 and result['ledger_error']<1e-10)
    return result,rows


def run(source, audit_source):
    start=time.monotonic()
    source=Path(source)
    saved=json.loads(source.with_name(source.name.replace('.dense.json','.path.json')).read_text())
    dense=json.loads(source.read_text())
    predecessor=json.loads(Path(audit_source).read_text())['result']
    if predecessor['source_sha256']!=hashlib.sha256(source.read_bytes()).hexdigest():
        raise ValueError('audit does not match dense source')
    worst_index=predecessor['worst']['segment']
    lo=max(0,worst_index-1);hi=min(len(dense['segments']),worst_index+2)
    old=[SimpleNamespace(t_old=s['t_old'],h=s['h'],y_old=np.asarray(s['y_old']),F=np.asarray(s['F'])) for s in dense['segments'][lo:hi]]
    problem=Problem(saved)
    a=old[0].t_old;b=old[-1].t_old+old[-1].h
    initial=old[0].y_old.copy()
    diagnostic=[]
    previous=None
    for k,seg in enumerate(old,lo):
        for t in np.linspace(seg.t_old,seg.t_old+seg.h,65):
            y=segment_value(seg,t);z=segment_derivative(seg,t)[:problem.p.N]
            residual=problem.residual(seg,t)
            exact_z,exact_f=problem.constrained(t,y[:problem.p.N])
            faces,speeds=problem.branches(t,y[:problem.p.N],z)
            topology=sorted((i,j) for i,j,_ in faces)
            switches=[]
            if previous is not None:
                switches=[key for key,v in speeds.items() if key in previous['speeds'] and v*previous['speeds'][key]<0]
            diagnostic.append(dict(segment=k,s=float(t),residual=residual,
                derivative_error=float(np.max(abs(z-exact_z))),
                solved_residual=float(np.max(abs(exact_f[:-4].reshape(problem.p.N,4)[:,0])/problem.p.mass)),
                faces=faces,speeds=speeds,sign_switches=switches,
                topology_changed=bool(previous is not None and topology!=previous['topology'])))
            previous=dict(speeds=speeds,topology=topology)
    Path('diagnostic.json').write_text(json.dumps(dict(window=[a,b],original_indices=[lo,hi],samples=diagnostic),indent=2))
    print(json.dumps(dict(stage='diagnostic_completed',window=[a,b],samples=len(diagnostic))),flush=True)
    for label,target in [('local_control',None),('residual_target_1e-9',1e-9)]:
        def progress(row):
            with open(label+'.attempts.jsonl','a') as f:f.write(json.dumps(row)+'\n')
            print(json.dumps(dict(case=label,**row)),flush=True)
        segments,end,attempts=integrate(problem.fun,initial,a,b,problem.check,target,progress=progress)
        Path(label+'.dense.json').write_text(json.dumps(dict(segments=[serialize(s) for s in segments],
            ts=[s.t_old for s in segments]+[b],window=[a,b],initial=initial.tolist(),end=end.tolist())))
        result,samples=audit_window(problem,segments)
        result.update(case=label,target=target,window=[a,b],attempts=len(attempts),
            rejections=sum(not q['accepted'] for q in attempts),
            endpoint_weight_difference=float(np.max(abs(end[:problem.p.N]-segment_value(old[-1],b)[:problem.p.N]))),
            local_only=True,spliced_into_global_path=False)
        Path(label+'.audit.json').write_text(json.dumps(dict(result=result,samples=samples),indent=2))
        with open('results.jsonl','a') as f:f.write(json.dumps(result)+'\n')
        print(json.dumps(result),flush=True)
    Path('completed.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,calls=problem.calls)))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('source');parser.add_argument('audit')
    args=parser.parse_args()
    run(args.source,args.audit)
