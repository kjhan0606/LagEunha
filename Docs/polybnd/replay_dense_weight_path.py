"""Audit saved DOP853 paths without integrating the weight ODE again.

Rebuild only the frozen predictor from the original source snapshot. Inspect
each segment separately, so both one-sided derivatives at knots are tested.
Rates are per normalized path time s, not physical time.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial import Polynomial as P
from dense_path_derivative import segment_derivative
from predicted_spacetime_flux import FluxAudit


def segment_value(segment, t):
    x = (t - segment.t_old) / segment.h
    value = np.zeros_like(segment.y_old)
    for i, f in enumerate(reversed(segment.F)):
        value = (value + f) * (x if i % 2 == 0 else 1 - x)
    return value + segment.y_old


def run(dense_file, output, nodes_per_segment=17, orders=(8, 16)):
    if nodes_per_segment < 3 or nodes_per_segment % 2 != 1:
        raise ValueError("nodes must be odd and >=3")
    low, high = orders
    if not 0 < low < high:
        raise ValueError("quadrature orders must increase")
    started = time.monotonic()
    dense_file = Path(dense_file)
    path_file = dense_file.with_name(dense_file.name.replace('.dense.json', '.path.json'))
    saved = json.loads(path_file.read_text())
    data = json.loads(dense_file.read_text())
    old = saved['result']
    audit = FluxAudit(old['kind'], old['n'], old['duration'])
    p = audit.path
    # Refuse a different predictor/generator implementation.
    generator_rebuild_difference = 0.
    for actual, expected in zip(p.x, saved['generator_polynomials']):
        for polynomial, coefficients in zip(actual, expected):
            generator_rebuild_difference = max(generator_rebuild_difference,
                float(np.max(abs(polynomial.coef-np.asarray(coefficients)))))
            np.testing.assert_allclose(polynomial.coef, coefficients, rtol=0, atol=1e-14)
    # Hardware may round predictor construction differently. Replay the original
    # generator exactly, and still enforce the independent old-sample check.
    p.x = [[P(coefficients) for coefficients in row] for row in saved['generator_polynomials']]
    segments = [SimpleNamespace(t_old=s['t_old'], h=s['h'],
                                y_old=np.asarray(s['y_old']), F=np.asarray(s['F']))
                for s in data['segments']]
    if len(segments) != len(data['ts']) - 1:
        raise ValueError('segment count mismatch')
    samples = []
    integral = {order: np.zeros(4*p.N+4) for order in orders}
    max_prefix = {order: 0. for order in integral}
    max_value_jump = 0.
    max_derivative_jump = 0.

    def flux(seg, t):
        y = segment_value(seg, t)
        p.w = [P([float(v)]) for v in y[:p.N]]
        z = segment_derivative(seg, t)[:p.N]
        result = audit.rate(t, z)
        if not np.all(np.isfinite(result)):
            raise ValueError('nonfinite flux')
        return result

    for k, seg in enumerate(segments):
        a, b = data['ts'][k:k+2]
        np.testing.assert_allclose([seg.t_old, seg.t_old+seg.h], [a, b], rtol=0, atol=2e-15)
        np.testing.assert_allclose(seg.y_old[:p.N], saved['accepted_weights'][k], rtol=0, atol=2e-15)
        if k:
            prior = segments[k-1]
            max_value_jump = max(max_value_jump, float(np.max(abs(segment_value(prior, a)-seg.y_old))))
            max_derivative_jump = max(max_derivative_jump, float(np.max(abs(
                segment_derivative(prior, a)[:p.N]-segment_derivative(seg, a)[:p.N]))))
        # Nested coarse/fine grids, including both sides of every knot.
        for j, t in enumerate(np.linspace(a, b, nodes_per_segment)):
            f = flux(seg, t)
            samples.append(dict(segment=k, node=j, s=float(t),
                                relative_mass_rate=float(np.max(abs(f[:-4].reshape(p.N,4)[:,0])/p.mass))))
        for order in integral:
            nodes, weights = np.polynomial.legendre.leggauss(order)
            for node, weight in zip(nodes, weights):
                integral[order] += seg.h * weight/2 * flux(seg, (a+b)/2+seg.h*node/2)
            max_prefix[order] = max(max_prefix[order], float(np.max(abs(
                integral[order][:-4].reshape(p.N,4)[:,0])/p.mass)))
        print(json.dumps(dict(segment=k+1, total=len(segments), source=str(dense_file))), flush=True)
    # Reproduce the old 21 samples before accepting denser audit results.
    reproduction_error = 0.
    for sample in saved['samples']:
        t = sample['t']
        k = int(np.clip(np.searchsorted(data['ts'], t, side='left')-1, 0, len(segments)-1))
        value = float(np.max(abs(flux(segments[k], t)[:-4].reshape(p.N,4)[:,0])/p.mass))
        reproduction_error = max(reproduction_error, abs(value-sample['dense_polynomial_relative_rate']))
    worst = max(samples, key=lambda s: s['relative_mass_rate'])
    f = integral[high][:-4].reshape(p.N, 4)
    delta_mass = (integral[high]-integral[low])[:-4].reshape(p.N,4)[:,0]
    result = dict(status='completed', source=str(dense_file),
                  source_sha256=hashlib.sha256(dense_file.read_bytes()).hexdigest(),
                  generator_rebuild_difference=generator_rebuild_difference,
                  seconds=time.monotonic()-started, segments=len(segments), samples=len(samples),
                  nodes_per_segment=nodes_per_segment, quadrature_orders=orders,
                  max_rate_coarse=max(s['relative_mass_rate'] for s in samples if s['node']%2==0),
                  max_rate_fine=worst['relative_mass_rate'], worst=worst,
                  integrated_relative_mass_defect=float(np.max(abs(f[:,0])/p.mass)),
                  quadrature_mass_difference=float(np.max(abs(delta_mass)/p.mass)),
                  quadrature_all_flux_difference=float(np.max(abs(integral[high]-integral[low]))),
                  max_prefix_relative_mass_defect=max_prefix,
                  boundary_ledger_error=float(np.max(abs(f.sum(0)-integral[high][-4:]))),
                  max_value_jump=max_value_jump, max_weight_derivative_jump=max_derivative_jump,
                  old_sample_reproduction_error=reproduction_error,
                  rate_units='per normalized time s; divide by duration for physical time',
                  duration=old['duration'], limitation='sampled check, not an all-time bound')
    result['gate'] = bool(result['max_rate_fine'] < 1e-8 and
                          result['integrated_relative_mass_defect'] < 1e-10 and
                          result['quadrature_mass_difference'] < 1e-11 and
                          result['boundary_ledger_error'] < 1e-10 and
                          reproduction_error < 1e-13 and max_value_jump < 1e-12)
    Path(output).write_text(json.dumps(dict(result=result, samples=samples,
        integrated_flux={str(k):v.tolist() for k,v in integral.items()}), indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('sources', nargs='+')
    parser.add_argument('--nodes', type=int, default=17)
    parser.add_argument('--orders', type=int, nargs=2, default=(8,16))
    args = parser.parse_args()
    for source in args.sources:
        case = Path(source).parent.name
        result = run(source, case+'.audit.json', args.nodes, tuple(args.orders))
        with open('results.jsonl', 'a') as stream:
            stream.write(json.dumps(dict(case=case, **result))+'\n')
        print(json.dumps(dict(case=case, **result)), flush=True)
