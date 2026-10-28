"""Volume-independent mass diagnostic with finite-width radial windows."""
import argparse
import json
from pathlib import Path
import numpy as np
from audit_kepler_native_radial_20261026 import ordered
from audit_kepler_t3_support_20261029 import digest


def analyze(initial, endpoint):
    a, m0 = ordered(initial)
    b, mt = ordered(endpoint)
    if m0['t'] != 0. or mt['t'] <= 0. or not np.array_equal(a['mass'], b['mass']):
        raise ValueError('Expected matching fixed masses and t0 reference')
    r0 = np.hypot(a['x']-3., a['y']-3.)
    rt = np.hypot(b['x']-3., b['y']-3.)
    if not np.isfinite(r0).all() or not np.isfinite(rt).all():
        raise ValueError('Nonfinite radius')
    foreign = r0 >= 2.75
    rows = []
    for lo, hi in [(.5,2.), (1.,2.), (1.5,2.)]:
        width = .02
        def smooth(r):
            return .25*(1.+np.tanh((r-lo)/width))*(1.-np.tanh((r-hi)/width))
        w0, wt = smooth(r0), smooth(rt)
        mass0 = np.sum(a['mass']*w0, dtype=np.longdouble)
        mass = np.sum(b['mass']*wt, dtype=np.longdouble)
        foreign0 = np.sum(a['mass'][foreign]*w0[foreign], dtype=np.longdouble)
        foreignt = np.sum(b['mass'][foreign]*wt[foreign], dtype=np.longdouble)
        rows.append({'rlo':lo, 'rhi':hi, 'smooth_width':width,
                     'initial_mass':float(mass0), 'current_mass':float(mass),
                     'change':float(mass-mass0), 'foreign_mass_change':float(foreignt-foreign0),
                     'nonforeign_mass_change':float(mass-mass0-foreignt+foreign0)})
    return {'initial':str(initial), 'endpoint':str(endpoint), 'time':mt['t'],
            'initial_sha256':digest(initial), 'endpoint_sha256':digest(endpoint), 'rows':rows,
            'scope':'Volume-independent smooth annulus mass using all generators and fixed masses. '
                    'Tanh width .02 avoids sensitivity of sharp cuts to polar rings exactly on r=.5 or2. '
                    'Not a face-integrated flux or claim that cell density errors equal this mass change.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('initial',type=Path)
    parser.add_argument('endpoint',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = analyze(args.initial,args.endpoint)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(result,indent=2))
