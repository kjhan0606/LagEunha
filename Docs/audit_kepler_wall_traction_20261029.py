"""Independent actual-wall traction quadrature versus native force-sum ledger."""
import argparse
import csv
import json
from pathlib import Path
import re
import numpy as np
from audit_kepler_restart_bytes_20261029 import compare


def analyze(root):
    config=json.loads((root/'stage.json').read_text());reference=Path(config['reference'])
    line=[s for s in (reference/'log').read_text().splitlines() if '[WALL_ACCEPT]' in s]
    if len(line)!=1:raise ValueError('Expected one accepted no-impact reference step')
    h=float(re.search(r'\bh=(\S+)',line[0])[1])
    match=re.search(r'wall_pressure_J=\(([^,]+),([^,]+),([^\)]+)\)',line[0])
    ledger=np.array([float(match[k]) for k in [1,2,3]])
    independent=np.zeros(3,dtype=np.longdouble);absolute=np.zeros(3,dtype=np.longdouble)
    stages=[];parity=[]
    for call,weight in enumerate(config['RK4_weights'],1):
        r=root/f'call{call}'
        if (r/'exit_status.txt').read_text().strip()!='0':raise ValueError('Native capture failed')
        with (r/f'face_dump_r0_call{call}.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        sides=set();force=np.zeros(3,dtype=np.longdouble);scale=0.;perimeter=0.;work_error=0.
        for q in rows:
            if int(q['call'])!=call or int(q['wall'])!=2:raise ValueError('Non-wall capture row')
            side=-1000-int(q['id_j']);sides.add(side)
            if not 0<=side<1024:raise ValueError('Invalid wall facet')
            pressure=float(q['pstar']);dx=float(q['dSx']);dy=float(q['dSy'])
            x=.5*(float(q['x1'])+float(q['x2']))-3.
            y=.5*(float(q['y1'])+float(q['y2']))-3.
            force+=np.array([pressure*dx,pressure*dy,pressure*(x*dy-y*dx)],dtype=np.longdouble)
            length=float(np.hypot(dx,dy));perimeter+=length;scale+=abs(pressure)*length
            work_error=max(work_error,abs(float(q['work'])+pressure*float(q['udotds'])))
        if len(sides)!=1024:raise ValueError('Incomplete actual wall coverage')
        independent+=h*weight*force;absolute+=h*weight*scale*np.array([1.,1.,2.75])
        stages.append({'call':call,'rows':len(rows),'wall_facets':len(sides),'wall_force_and_torque':force.astype(float).tolist(),
                       'wall_perimeter':perimeter,'maximum_relative_sweep_work_formula_error':work_error})
        parity.append(compare(reference/'khout.011581.dat',r/'khout.011581.dat'))
    delta=np.abs(independent-ledger)
    tolerance=config['relative_tolerance_on_absolute_wall_impulse_scale']*absolute
    checks={'actual_wall_traction_matches_force_sum_ledger':bool(np.all(delta<=tolerance)),
            'all_captures_noninvasive_full_restart_exact':all(p['passed'] for p in parity),
            'wall_work_is_geometric_relative_sweep_correction':all(s['maximum_relative_sweep_work_formula_error']<1.e-18 for s in stages)}
    return {'passed':all(checks.values()),'checks':checks,'h':h,'stages':stages,
            'wall_impulse_independent_faces':independent.astype(float).tolist(),'wall_impulse_native_ledger':ledger.tolist(),
            'absolute_difference':delta.astype(float).tolist(),'declared_tolerance':tolerance.astype(float).tolist(),
            'capture_full_restart_comparisons':parity,
            'scope':'All actual boundary faces at four actual native RK4 stages of the same no-impact step. Separately evaluates integral P n dA and torque; not derived from fluid acceleration. The wall is stationary: work column is the relative-sweep energy correction, not external wall power. No new physical evolution or long-time claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=analyze(args.root)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r,indent=2))
    if not r['passed']:raise SystemExit(1)
