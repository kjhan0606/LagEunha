"""Independent fixed-step refinements of the completed Q extension cases."""
import argparse
import json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from coupled_factor_path_audit import run

CASES = {'wave_n4': (4, .01), 'wave_strong_n3': (3, .1)}
COUNTS = (16, 32, 64)
GROUPS = {'sites': slice(0, 2), 'weights': slice(2, 3),
          'q': slice(3, 5), 'Q': slice(5, 9), 'p': slice(9, 11),
          'D': slice(11, 15), 'Pi': slice(15, 19), 'E': slice(19, 20)}


def compare(states, reference):
    result = {}
    for name, sl in GROUPS.items():
        differences = [float(np.linalg.norm(a[:, sl]-b[:, sl]))
                       for a, b in zip(states[:-1], states[1:])]
        errors = [float(np.linalg.norm(a[:, sl]-reference[:, sl])) for a in states]
        resolved = min(differences) > 1e-12
        order = float(np.log2(differences[0]/differences[1])) if resolved else None
        result[name] = dict(successive_L2_differences=differences,
            adaptive_reference_L2_errors=errors, observed_order=order,
            resolved=resolved, second_order_gate=(1.5 <= order <= 2.5) if resolved else None)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--steps', type=int, choices=COUNTS)
    parser.add_argument('--report', action='store_true')
    args = parser.parse_args()
    predecessor = json.loads(Path('predecessor_completion.json').read_text())
    if not predecessor.get('passed'):
        raise RuntimeError('Q predecessor failed')
    if not args.report:
        if args.case is None or args.steps is None:
            parser.error('--case and --steps required')
        n, amplitude = CASES[args.case]
        if not predecessor['cases'][args.case]['passed']:
            raise RuntimeError('case predecessor failed')
        # Existing directories are deliberately refused: never restart over saved data.
        result = run(args.steps, n=n, amplitude=amplitude)
        atomic('case_summary.json', dict(case=args.case, amplitude=amplitude, result=result))
        raise SystemExit(0 if result.get('state_gate') else 2)
    summaries = {}
    for name, (n, amplitude) in CASES.items():
        try:
            rows = [json.loads((Path(name)/str(k)/'case_summary.json').read_text()) for k in COUNTS]
            if not all(r['case'] == name and r['amplitude'] == amplitude and
                       r['result'].get('state_gate') for r in rows):
                raise RuntimeError('missing/failed state gate or case mismatch')
            checkpoints = [json.loads((Path(name)/str(k)/f'steps_{k}'/'checkpoint.json').read_text()) for k in COUNTS]
            reference = json.loads((Path('references')/f'{name}.json').read_text())
            cfg = reference['config']
            if (cfg['n'], cfg['amplitude'], cfg['end'], reference['t']) != (n, amplitude, .05, .05):
                raise RuntimeError('reference configuration/time mismatch')
            if not all(c['time'] == .05 and c['step'] == k for c, k in zip(checkpoints, COUNTS)):
                raise RuntimeError('incomplete endpoint')
            for c in checkpoints:
                np.testing.assert_allclose(c['mass'], reference['mass'], rtol=0, atol=1e-15)
            states = [np.asarray(c['state'])[:-4].reshape(n*n, 20) for c in checkpoints]
            ref = np.asarray(reference['state'])[:-4].reshape(n*n, 20)
            groups = compare(states, ref)
            resolved = [v for v in groups.values() if v['resolved']]
            passed = bool(resolved) and all(v['second_order_gate'] for v in resolved)
            summaries[name] = dict(passed=passed, groups=groups,
                fixed_step_results=[r['result'] for r in rows])
        except Exception as exc:
            summaries[name] = dict(passed=False, error=repr(exc))
    report = dict(cases=summaries, passed=all(v['passed'] for v in summaries.values()),
        note='Temporal self-convergence at fixed spatial resolution. Adaptive solution is a numerical comparison, not an exact solution. Coarse-step mass path gates are reported separately.')
    atomic('completion_report.json', report)
    print(json.dumps(report), flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__ == '__main__':
    main()
